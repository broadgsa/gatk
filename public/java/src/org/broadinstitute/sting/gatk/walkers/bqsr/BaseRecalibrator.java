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

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.classloader.GATKLiteUtils;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.recalibration.*;
import org.broadinstitute.sting.utils.recalibration.covariates.Covariate;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * First pass of the base quality score recalibration -- Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
 *
 * <p>
 * This walker is designed to work as the first pass in a two-pass processing step. It does a by-locus traversal operating
 * only at sites that are not in dbSNP. We assume that all reference mismatches we see are therefore errors and indicative
 * of poor base quality. This walker generates tables based on various user-specified covariates (such as read group,
 * reported quality score, cycle, and context). Since there is a large amount of data one can then calculate an empirical
 * probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations.
 * The output file is a table (of the several covariate values, num observations, num mismatches, empirical quality score).
 * <p>
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and will be added for the user regardless of whether or not they were specified.
 *
 * <p>
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
 *   -T BaseRecalibrator \
 *   -I my_reads.bam \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -o recal_data.grp
 * </pre>
 */

@DocumentedGATKFeature(groupName = "BAM Processing and Analysis Tools", extraDocs = {CommandLineGATK.class})
@BAQMode(ApplicationTime = ReadTransformer.ApplicationTime.FORBIDDEN)
@ReadFilters({MappingQualityZeroFilter.class, MappingQualityUnavailableFilter.class, UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class})
@PartitionBy(PartitionType.READ)
public class BaseRecalibrator extends ReadWalker<Long, Long> implements NanoSchedulable {
    @ArgumentCollection
    private final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection(); // all the command line arguments for BQSR and it's covariates

    @Advanced
    @Argument(fullName = "bqsrBAQGapOpenPenalty", shortName="bqsrBAQGOP", doc="BQSR BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets", required = false)
    public double BAQGOP = BAQ.DEFAULT_GOP;

    private QuantizationInfo quantizationInfo; // an object that keeps track of the information necessary for quality score quantization

    private RecalibrationTables recalibrationTables;

    private Covariate[] requestedCovariates; // list to hold the all the covariate objects that were requested (required + standard + experimental)

    private RecalibrationEngine recalibrationEngine;

    private int minimumQToUse;

    protected static final String COVARS_ATTRIBUTE = "COVARS"; // used to store covariates array as a temporary attribute inside GATKSAMRecord.\

    private static final String NO_DBSNP_EXCEPTION = "This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.";

    private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
    private IndexedFastaSequenceFile referenceReader; // fasta reference reader for use with BAQ calculation

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    public void initialize() {

        baq = new BAQ(BAQGOP); // setup the BAQ object with the provided gap open penalty

        // check for unsupported access
        if (getToolkit().isGATKLite() && !getToolkit().getArguments().disableIndelQuals)
            throw new UserException.NotSupportedInGATKLite("base insertion/deletion recalibration is not supported, please use the --disable_indel_quals argument");

        if (RAC.FORCE_PLATFORM != null)
            RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM;

        if (RAC.knownSites.isEmpty() && !RAC.RUN_WITHOUT_DBSNP) // Warn the user if no dbSNP file or other variant mask was specified
            throw new UserException.CommandLineException(NO_DBSNP_EXCEPTION);

        if (RAC.LIST_ONLY) {
            RecalUtils.listAvailableCovariates(logger);
            System.exit(0);
        }
        RAC.existingRecalibrationReport = getToolkit().getArguments().BQSR_RECAL_FILE; // if we have a recalibration file, record it so it goes on the report table

        Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = RecalUtils.initializeCovariates(RAC); // initialize the required and optional covariates
        ArrayList<Covariate> requiredCovariates = covariates.getFirst();
        ArrayList<Covariate> optionalCovariates = covariates.getSecond();

        requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
        int covariateIndex = 0;
        for (final Covariate covariate : requiredCovariates)
            requestedCovariates[covariateIndex++] = covariate;
        for (final Covariate covariate : optionalCovariates)
            requestedCovariates[covariateIndex++] = covariate;

        logger.info("The covariates being used here: ");
        for (Covariate cov : requestedCovariates) { // list all the covariates being used
            logger.info("\t" + cov.getClass().getSimpleName());
            cov.initialize(RAC); // initialize any covariate member variables using the shared argument collection
        }

        try {
            RAC.RECAL_TABLE = new PrintStream(RAC.RECAL_TABLE_FILE);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(RAC.RECAL_TABLE_FILE, e);
        }

        int numReadGroups = 0;
        for ( final SAMFileHeader header : getToolkit().getSAMFileHeaders() )
            numReadGroups += header.getReadGroups().size();
        recalibrationTables = new RecalibrationTables(requestedCovariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG);

        recalibrationEngine = initializeRecalibrationEngine();
        recalibrationEngine.initialize(requestedCovariates, recalibrationTables);

        minimumQToUse = getToolkit().getArguments().PRESERVE_QSCORES_LESS_THAN;

        try {
            // fasta reference reader for use with BAQ calculation
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile, e);
        }

    }

    private RecalibrationEngine initializeRecalibrationEngine() {

        final Class recalibrationEngineClass = GATKLiteUtils.getProtectedClassIfAvailable(RecalibrationEngine.class);
        try {
            final Constructor constructor = recalibrationEngineClass.getDeclaredConstructor((Class[])null);
            constructor.setAccessible(true);
            return (RecalibrationEngine)constructor.newInstance();
        }
        catch (Exception e) {
            throw new ReviewedStingException("Unable to create RecalibrationEngine class instance " + recalibrationEngineClass.getSimpleName());
        }
    }

    private boolean isLowQualityBase( final GATKSAMRecord read, final int offset ) {
        return read.getBaseQualities()[offset] < minimumQToUse;
    }

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     * whether or not the base matches the reference at this particular location
     */
    public Long map( final ReferenceContext ref, final GATKSAMRecord originalRead, final RefMetaDataTracker metaDataTracker ) {

        final GATKSAMRecord read = ReadClipper.hardClipSoftClippedBases( ReadClipper.hardClipAdaptorSequence(originalRead) );
        if( read.isEmpty() ) { return 0L; } // the whole read was inside the adaptor so skip it

        RecalUtils.parsePlatformForRead(read, RAC);
        if (!RecalUtils.isColorSpaceConsistent(RAC.SOLID_NOCALL_STRATEGY, read)) { // parse the solid color space and check for color no-calls
            return 0L; // skip this read completely
        }
        read.setTemporaryAttribute(COVARS_ATTRIBUTE, RecalUtils.computeCovariates(read, requestedCovariates));

        final boolean[] skip = calculateSkipArray(read, metaDataTracker); // skip known sites of variation as well as low quality and non-regular bases
        final int[] isSNP = calculateIsSNP(read, ref, originalRead);
        final int[] isInsertion = calculateIsIndel(read, EventType.BASE_INSERTION);
        final int[] isDeletion = calculateIsIndel(read, EventType.BASE_DELETION);
        final byte[] baqArray = calculateBAQArray(read);

        if( baqArray != null ) { // some reads just can't be BAQ'ed
            final double[] snpErrors = calculateFractionalErrorArray(isSNP, baqArray);
            final double[] insertionErrors = calculateFractionalErrorArray(isInsertion, baqArray);
            final double[] deletionErrors = calculateFractionalErrorArray(isDeletion, baqArray);
            recalibrationEngine.updateDataForRead(read, skip, snpErrors, insertionErrors, deletionErrors);
            return 1L;
        } else {
            return 0L;
        }
    }

    protected boolean[] calculateSkipArray( final GATKSAMRecord read, final RefMetaDataTracker metaDataTracker ) {
        final byte[] bases = read.getReadBases();
        final boolean[] skip = new boolean[bases.length];
        final boolean[] knownSites = calculateKnownSites(read, metaDataTracker.getValues(RAC.knownSites));
        for( int iii = 0; iii < bases.length; iii++ ) {
            skip[iii] = !BaseUtils.isRegularBase(bases[iii]) || isLowQualityBase(read, iii) || knownSites[iii] || badSolidOffset(read, iii);
        }
        return skip;
    }

    protected boolean badSolidOffset( final GATKSAMRecord read, final int offset ) {
        return ReadUtils.isSOLiDRead(read) && RAC.SOLID_RECAL_MODE != RecalUtils.SOLID_RECAL_MODE.DO_NOTHING && !RecalUtils.isColorSpaceConsistent(read, offset);
    }

    protected boolean[] calculateKnownSites( final GATKSAMRecord read, final List<Feature> features ) {
        final int readLength = read.getReadBases().length;
        final boolean[] knownSites = new boolean[readLength];
        Arrays.fill(knownSites, false);
        for( final Feature f : features ) {
            int featureStartOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), f.getStart(), ReadUtils.ClippingTail.LEFT_TAIL, true); // BUGBUG: should I use LEFT_TAIL here?
            if( featureStartOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureStartOnRead = 0;
            }

            int featureEndOnRead = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), f.getEnd(), ReadUtils.ClippingTail.LEFT_TAIL, true);
            if( featureEndOnRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
                featureEndOnRead = readLength;
            }

            if( featureStartOnRead > readLength ) {
                featureStartOnRead = featureEndOnRead = readLength;
            }

            Arrays.fill(knownSites, Math.max(0, featureStartOnRead), Math.min(readLength, featureEndOnRead + 1), true);
        }
        return knownSites;
    }

    // BUGBUG: can be merged with calculateIsIndel
    protected static int[] calculateIsSNP( final GATKSAMRecord read, final ReferenceContext ref, final GATKSAMRecord originalRead ) {
        final byte[] readBases = read.getReadBases();
        final byte[] refBases = Arrays.copyOfRange(ref.getBases(), read.getAlignmentStart() - originalRead.getAlignmentStart(), ref.getBases().length + read.getAlignmentEnd() - originalRead.getAlignmentEnd());
        final int[] snp = new int[readBases.length];
        int readPos = 0;
        int refPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        snp[readPos] = ( BaseUtils.basesAreEqual(readBases[readPos], refBases[refPos]) ? 0 : 1 );
                        readPos++;
                        refPos++;
                    }
                    break;
                case D:
                case N:
                    refPos += elementLength;
                    break;
                case I:
                case S: // ReferenceContext doesn't have the soft clipped bases!
                    readPos += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new ReviewedStingException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return snp;
    }

    protected static int[] calculateIsIndel( final GATKSAMRecord read, final EventType mode ) {
        final byte[] readBases = read.getReadBases();
        final int[] indel = new int[readBases.length];
        Arrays.fill(indel, 0);
        int readPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                case S:
                {
                    readPos += elementLength;
                    break;
                }
                case D:
                {
                    final int index = ( read.getReadNegativeStrandFlag() ? readPos : ( readPos > 0 ? readPos - 1 : readPos ) );
                    indel[index] = ( mode.equals(EventType.BASE_DELETION) ? 1 : 0 );
                    break;
                }
                case I:
                {
                    final boolean forwardStrandRead = !read.getReadNegativeStrandFlag();
                    if( forwardStrandRead ) {
                        indel[(readPos > 0 ? readPos - 1 : readPos)] = ( mode.equals(EventType.BASE_INSERTION) ? 1 : 0 );
                    }
                    for (int iii = 0; iii < elementLength; iii++) {
                        readPos++;
                    }
                    if( !forwardStrandRead ) {
                        indel[(readPos < indel.length ? readPos : readPos - 1)] = ( mode.equals(EventType.BASE_INSERTION) ? 1 : 0 );
                    }
                    break;
                }
                case N:
                case H:
                case P:
                    break;
                default:
                    throw new ReviewedStingException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return indel;
    }

    protected static double[] calculateFractionalErrorArray( final int[] errorArray, final byte[] baqArray ) {
        if(errorArray.length != baqArray.length ) {
            throw new ReviewedStingException("Array length mismatch detected. Malformed read?");
        }

        final byte NO_BAQ_UNCERTAINTY = (byte)'@';
        final int BLOCK_START_UNSET = -1;

        final double[] fractionalErrors = new double[baqArray.length];
        Arrays.fill(fractionalErrors, 0.0);
        boolean inBlock = false;
        int blockStartIndex = BLOCK_START_UNSET;
        int iii;
        for( iii = 0; iii < fractionalErrors.length; iii++ ) {
            if( baqArray[iii] == NO_BAQ_UNCERTAINTY ) {
                if( !inBlock ) {
                    fractionalErrors[iii] = (double) errorArray[iii];
                } else {
                    calculateAndStoreErrorsInBlock(iii, blockStartIndex, errorArray, fractionalErrors);
                    inBlock = false; // reset state variables
                    blockStartIndex = BLOCK_START_UNSET; // reset state variables
                }
            } else {
                inBlock = true;
                if( blockStartIndex == BLOCK_START_UNSET ) { blockStartIndex = iii; }
            }
        }
        if( inBlock ) {
            calculateAndStoreErrorsInBlock(iii-1, blockStartIndex, errorArray, fractionalErrors);
        }
        if( fractionalErrors.length != errorArray.length ) {
            throw new ReviewedStingException("Output array length mismatch detected. Malformed read?");
        }
        return fractionalErrors;
    }

    private static void calculateAndStoreErrorsInBlock( final int iii,
                                                        final int blockStartIndex,
                                                        final int[] errorArray,
                                                        final double[] fractionalErrors ) {
        int totalErrors = 0;
        for( int jjj = Math.max(0,blockStartIndex-1); jjj <= iii; jjj++ ) {
            totalErrors += errorArray[jjj];
        }
        for( int jjj = Math.max(0, blockStartIndex-1); jjj <= iii; jjj++ ) {
            fractionalErrors[jjj] = ((double) totalErrors) / ((double)(iii - Math.max(0,blockStartIndex-1) + 1));
        }
    }

    private byte[] calculateBAQArray( final GATKSAMRecord read ) {
        baq.baqRead(read, referenceReader, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
        return BAQ.getBAQTag(read);
    }

    /**
     * Initialize the reduce step by returning 0L
     *
     * @return returns 0L
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

    @Override
    public void onTraversalDone(Long result) {
        logger.info("Calculating quantized quality scores...");
        quantizeQualityScores();

        logger.info("Writing recalibration report...");
        generateReport();
        logger.info("...done!");

        if ( RAC.RECAL_PDF_FILE != null ) {
            logger.info("Generating recalibration plots...");
            generatePlots();
        }

        logger.info("Processed: " + result + " reads");
    }

    private void generatePlots() {
        File recalFile = getToolkit().getArguments().BQSR_RECAL_FILE;
        if (recalFile != null) {
            RecalibrationReport report = new RecalibrationReport(recalFile);
            RecalUtils.generateRecalibrationPlot(RAC, report.getRecalibrationTables(), recalibrationTables, requestedCovariates);
        }
        else
            RecalUtils.generateRecalibrationPlot(RAC, recalibrationTables, requestedCovariates);
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
        RecalUtils.outputRecalibrationReport(RAC, quantizationInfo, recalibrationTables, requestedCovariates);
    }
}