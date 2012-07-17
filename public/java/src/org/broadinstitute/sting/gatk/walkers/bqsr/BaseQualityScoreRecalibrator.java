/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.bqsr;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.classloader.ProtectedPackageSource;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.recalibration.RecalibrationTables;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.*;

/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various user-specified covariates (such as reported quality score, cycle, and dinucleotide).
 *
 * <p>
 * This walker is designed to work as the first pass in a two-pass processing step. It does a by-locus traversal operating
 * only at sites that are not in dbSNP. We assume that all reference mismatches we see are therefore errors and indicative
 * of poor base quality. This walker generates tables based on various user-specified covariates (such as read group,
 * reported quality score, cycle, and dinucleotide). Since there is a large amount of data one can then calculate an empirical
 * probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.
 * The output file is a CSV list of (the several covariate values, num observations, num mismatches, empirical quality score).
 * <p>
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and will be added for the user regardless of whether or not they were specified.
 *
 * <p>
 * See the GATK wiki for a tutorial and example recalibration accuracy plots.
 * http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
 *
 * <h2>Input</h2>
 * <p>
 * The input read data whose base quality scores need to be assessed.
 * <p>
 * A database of known polymorphic sites to skip over.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A GATK Report file with many tables:
 * <ol>
 *     <li>The list of arguments</li>
 *     <li>The quantized qualities table</li>
 *     <li>The recalibration table by read group</li>
 *     <li>The recalibration table by quality score</li>
 *     <li>The recalibration table for all the optional covariates</li>
 * </ol>
 *
 * The GATK Report is intended to be easy to read by humans or computers. Check out the documentation of the GATKReport to learn how to manipulate this table.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -T BaseQualityScoreRecalibrator \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -o recal_data.grp
 * </pre>
 */

@BAQMode(ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
@By(DataSource.READS)

@ReadFilters({MappingQualityZeroFilter.class, MappingQualityUnavailableFilter.class})                                   // only look at covered loci, not every loci of the reference file
@Requires({DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES})                                         // filter out all reads with zero or unavailable mapping quality
@PartitionBy(PartitionType.LOCUS)                                                                                       // this walker requires both -I input.bam and -R reference.fasta

public class BaseQualityScoreRecalibrator extends LocusWalker<Long, Long> implements TreeReducible<Long> {
    @ArgumentCollection
    private final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();                          // all the command line arguments for BQSR and it's covariates

    private QuantizationInfo quantizationInfo;                                                                          // an object that keeps track of the information necessary for quality score quantization 
    
    private RecalibrationTables recalibrationTables;
    
    private Covariate[] requestedCovariates;                                                                            // list to hold the all the covariate objects that were requested (required + standard + experimental)

    private RecalibrationEngine recalibrationEngine;

    private int minimumQToUse;

    protected static final String SKIP_RECORD_ATTRIBUTE = "SKIP";                                                         // used to label reads that should be skipped.
    protected static final String SEEN_ATTRIBUTE = "SEEN";                                                                // used to label reads as processed.
    protected static final String COVARS_ATTRIBUTE = "COVARS";                                                            // used to store covariates array as a temporary attribute inside GATKSAMRecord.\

    private static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";


    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    public void initialize() {

        // check for unsupported access
        if (getToolkit().isGATKLite() && !getToolkit().getArguments().disableIndelQuals)
            throw new UserException.NotSupportedInGATKLite("base insertion/deletion recalibration is not supported, please use the --disable_indel_quals argument");

        if (RAC.FORCE_PLATFORM != null)
            RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;

        if (RAC.knownSites.isEmpty() && !RAC.RUN_WITHOUT_DBSNP)                                                         // Warn the user if no dbSNP file or other variant mask was specified
            throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);

        if (RAC.LIST_ONLY) {
            RecalDataManager.listAvailableCovariates(logger);
            System.exit(0);
        }
        RAC.recalibrationReport = getToolkit().getArguments().BQSR_RECAL_FILE;                                          // if we have a recalibration file, record it so it goes on the report table

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalDataManager.initializeCovariates(RAC);       // initialize the required and optional covariates
        ArrayList<Covariate> requiredCovariates = covariates.getFirst();
        ArrayList<Covariate> optionalCovariates = covariates.getSecond();

        requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        int covariateIndex = 0;
        for (final Covariate covariate : requiredCovariates)
            requestedCovariates[covariateIndex++] = covariate;
        for (final Covariate covariate : optionalCovariates)
            requestedCovariates[covariateIndex++] = covariate;

        logger.info("The covariates being used here: ");
        for (Covariate cov : requestedCovariates) {                                                                     // list all the covariates being used
            logger.info("\t" + cov.getClass().getSimpleName());
            cov.initialize(RAC);                                                                                        // initialize any covariate member variables using the shared argument collection
        }

        int numReadGroups = 0;
        for ( final SAMFileHeader header : getToolkit().getSAMFileHeaders() )
            numReadGroups += header.getReadGroups().size();
        recalibrationTables = new RecalibrationTables(requestedCovariates, numReadGroups);

        recalibrationEngine = initializeRecalibrationEngine();
        recalibrationEngine.initialize(requestedCovariates, recalibrationTables);

        minimumQToUse = getToolkit().getArguments().PRESERVE_QSCORES_LESS_THAN;
    }

    private RecalibrationEngine initializeRecalibrationEngine() {
        List<Class<? extends RecalibrationEngine>> REclasses = new PluginManager<RecalibrationEngine>(RecalibrationEngine.class).getPlugins();
        if ( REclasses.isEmpty() )
            throw new ReviewedStingException("The RecalibrationEngine class is not found; repository must be corrupted");

        Class c = null;
        for ( Class<? extends RecalibrationEngine> REclass : REclasses ) {
            if ( REclass.isAssignableFrom(ProtectedPackageSource.class) ) {
                c = REclass;
                break;
            }
        }
        if ( c == null )
            c = REclasses.get(0);

        try {
            Constructor constructor = c.getDeclaredConstructor((Class[])null);
            constructor.setAccessible(true);
            return (RecalibrationEngine)constructor.newInstance();
        }
        catch (Exception e) {
            throw new ReviewedStingException("Unable to create RecalibrationEngine class instance " + c.getSimpleName());
        }
    }

    private boolean readHasBeenSkipped(GATKSAMRecord read) {
        return read.containsTemporaryAttribute(SKIP_RECORD_ATTRIBUTE);
    }

    private boolean isLowQualityBase(GATKSAMRecord read, int offset) {
        return read.getBaseQualities()[offset] < minimumQToUse;
    }

    private boolean readNotSeen(GATKSAMRecord read) {
        return !read.containsTemporaryAttribute(SEEN_ATTRIBUTE);
    }

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     *
     * @param tracker the reference metadata tracker
     * @param ref     the reference context
     * @param context the alignment context
     * @return returns 1, but this value isn't used in the reduce step
     */
    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        long countedSites = 0L;
        if (tracker.getValues(RAC.knownSites).size() == 0) {                                                            // Only analyze sites not present in the provided known sites
            for (final PileupElement p : context.getBasePileup()) {
                final GATKSAMRecord read = p.getRead();
                final int offset = p.getOffset();

                if (readHasBeenSkipped(read) || isLowQualityBase(read, offset))                                         // This read has been marked to be skipped or base is low quality (we don't recalibrate low quality bases)
                    continue;

                if (readNotSeen(read)) {
                    read.setTemporaryAttribute(SEEN_ATTRIBUTE, true);
                    RecalDataManager.parsePlatformForRead(read, RAC);
                    if (RecalDataManager.isColorSpaceConsistent(RAC.SOLID_NOCALL_STRATEGY, read)) {
                        read.setTemporaryAttribute(SKIP_RECORD_ATTRIBUTE, true);
                        continue;
                    }
                    read.setTemporaryAttribute(COVARS_ATTRIBUTE, RecalDataManager.computeCovariates(read, requestedCovariates));
                }

                if (!ReadUtils.isSOLiDRead(read) ||                                                                     // SOLID bams have inserted the reference base into the read if the color space in inconsistent with the read base so skip it
                    RAC.SOLID_RECAL_MODE == RecalDataManager.SOLID_RECAL_MODE.DO_NOTHING ||
                        RecalDataManager.isColorSpaceConsistent(read, offset))
                    recalibrationEngine.updateDataForPileupElement(p, ref.getBase());                                                             // This base finally passed all the checks for a good base, so add it to the big data hashmap
            }
            countedSites++;
        }

        return countedSites;
    }

    /**
     * Initialize the reduce step by creating a PrintStream from the filename specified as an argument to the walker.
     *
     * @return returns A PrintStream created from the -recalFile filename argument specified to the walker
     */
    public Long reduceInit() {
        return 0L;
    }

    /**
     * The Reduce method doesn't do anything for this walker.
     *
     * @param mapped Result of the map. This value is immediately ignored.
     * @param sum    The summing CountedData used to output the CSV data
     * @return returns The sum used to output the CSV data
     */
    public Long reduce(Long mapped, Long sum) {
        sum += mapped;
        return sum;
    }

    public Long treeReduce(Long sum1, Long sum2) {
        sum1 += sum2;
        return sum1;
    }

    @Override
    public void onTraversalDone(Long result) {
        logger.info("Calculating quantized quality scores...");
        quantizeQualityScores();
        if (!RAC.NO_PLOTS) {
            logger.info("Generating recalibration plots...");
            generatePlots();
        }
        logger.info("Writing recalibration report...");
        generateReport();
        logger.info("...done!");
        logger.info("Processed: " + result + " sites");
    }

    private void generatePlots() {
        File recalFile = getToolkit().getArguments().BQSR_RECAL_FILE;
        if (recalFile != null) {
            RecalibrationReport report = new RecalibrationReport(recalFile);
            RecalDataManager.generateRecalibrationPlot(RAC.RECAL_FILE, report.getRecalibrationTables(), recalibrationTables, requestedCovariates, RAC.KEEP_INTERMEDIATE_FILES);
        }
        else
            RecalDataManager.generateRecalibrationPlot(RAC.RECAL_FILE, recalibrationTables, requestedCovariates, RAC.KEEP_INTERMEDIATE_FILES);
    }


    /**
     * go through the quality score table and use the # observations and the empirical quality score
     * to build a quality score histogram for quantization. Then use the QuantizeQual algorithm to
     * generate a quantization map (recalibrated_qual -> quantized_qual)
     */
    private void quantizeQualityScores() {
        quantizationInfo = new QuantizationInfo(recalibrationTables, RAC.QUANTIZING_LEVELS);
    }

    private void generateReport() {
        PrintStream output;
        try {
            output = new PrintStream(RAC.RECAL_FILE);
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(RAC.RECAL_FILE, "could not be created");
        }

        RecalDataManager.outputRecalibrationReport(RAC, quantizationInfo, recalibrationTables, requestedCovariates, output);
    }
}

