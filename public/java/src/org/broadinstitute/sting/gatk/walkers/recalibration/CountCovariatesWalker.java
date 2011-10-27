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

package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

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
 * A recalibration table file in CSV format that is used by the TableRecalibration walker.
 * It is a comma-separated text file relating the desired covariates to the number of such bases and their rate of mismatch in the genome, and its implied empirical quality score.  
 *
 * The first 20 lines of such a file is shown below.  
 * * The file begins with a series of comment lines describing:
 * ** The number of counted loci
 * ** The number of counted bases
 * ** The number of skipped loci and the fraction skipped, due to presence in dbSNP or bad reference bases
 * 
 * * After the comments appears a header line indicating which covariates were used as well as the ordering of elements in the subsequent records.  
 *
 * * After the header, data records occur one per line until the end of the file. The first several items on a line are the values of the individual covariates and will change
 * depending on which covariates were specified at runtime. The last three items are the data- that is, number of observations for this combination of covariates, number of 
 * reference mismatches, and the raw empirical quality score calculated by phred-scaling the mismatch rate.
 * 
 * <pre>
 * # Counted Sites    19451059
 * # Counted Bases    56582018
 * # Skipped Sites    82666
 * # Fraction Skipped 1 / 235 bp
 * ReadGroup,QualityScore,Cycle,Dinuc,nObservations,nMismatches,Qempirical
 * SRR006446,11,65,CA,9,1,10
 * SRR006446,11,48,TA,10,0,40
 * SRR006446,11,67,AA,27,0,40
 * SRR006446,11,61,GA,11,1,10
 * SRR006446,12,34,CA,47,1,17
 * SRR006446,12,30,GA,52,1,17
 * SRR006446,12,36,AA,352,1,25
 * SRR006446,12,17,TA,182,11,12
 * SRR006446,11,48,TG,2,0,40
 * SRR006446,11,67,AG,1,0,40
 * SRR006446,12,34,CG,9,0,40
 * SRR006446,12,30,GG,43,0,40
 * ERR001876,4,31,AG,1,0,40
 * ERR001876,4,31,AT,2,2,1
 * ERR001876,4,31,CA,1,0,40
 * </pre>
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx4g -jar GenomeAnalysisTK.jar \
 *   -R resources/Homo_sapiens_assembly18.fasta \
 *   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
 *   -knownSites another/optional/setOfSitesToMask.vcf \
 *   -I my_reads.bam \
 *   -T CountCovariates \
 *   -cov ReadGroupCovariate \
 *   -cov QualityScoreCovariate \
 *   -cov CycleCovariate \
 *   -cov DinucCovariate \
 *   -recalFile my_reads.recal_data.csv
 * </pre>
 *
 */

@BAQMode(ApplicationTime = BAQ.ApplicationTime.FORBIDDEN)
@By( DataSource.READS ) // Only look at covered loci, not every loci of the reference file
@ReadFilters( {MappingQualityZeroFilter.class, MappingQualityUnavailableFilter.class} ) // Filter out all reads with zero or unavailable mapping quality
@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta
@PartitionBy(PartitionType.LOCUS)
public class CountCovariatesWalker extends LocusWalker<CountCovariatesWalker.CountedData, CountCovariatesWalker.CountedData> implements TreeReducible<CountCovariatesWalker.CountedData> {

    /////////////////////////////
    // Constants
    /////////////////////////////
    private static final String SKIP_RECORD_ATTRIBUTE = "SKIP"; //used to label GATKSAMRecords that should be skipped.
    private static final String SEEN_ATTRIBUTE = "SEEN"; //used to label GATKSAMRecords as processed.
    private static final String COVARS_ATTRIBUTE = "COVARS"; //used to store covariates array as a temporary attribute inside GATKSAMRecord.

    /////////////////////////////
    // Shared Arguments
    /////////////////////////////
    @ArgumentCollection private RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    /**
     * This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference,
     * so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites.
     */
    @Input(fullName="knownSites", shortName = "knownSites", doc="A database of known polymorphic sites to skip over in the recalibration algorithm", required=false)
    public List<RodBinding<Feature>> knownSites = Collections.emptyList();

    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.
     */
    @Output(fullName="recal_file", shortName="recalFile", required=true, doc="Filename for the output covariates table recalibration file")
    @Gather(CountCovariatesGatherer.class)
    public PrintStream RECAL_FILE;

    @Argument(fullName="list", shortName="ls", doc="List the available covariates and exit", required=false)
    private boolean LIST_ONLY = false;

    /**
     * See the -list argument to view available covariates.
     */
    @Argument(fullName="covariate", shortName="cov", doc="Covariates to be used in the recalibration. Each covariate is given as a separate cov parameter. ReadGroup and ReportedQuality are required covariates and are already added for you.", required=false)
    private String[] COVARIATES = null;
    @Argument(fullName="standard_covs", shortName="standard", doc="Use the standard set of covariates in addition to the ones listed using the -cov argument", required=false)
    private boolean USE_STANDARD_COVARIATES = false;

    /////////////////////////////
    // Debugging-only Arguments
    /////////////////////////////
    @Argument(fullName="dont_sort_output", shortName="unsorted", required=false, doc="If specified, the output table recalibration csv file will be in an unsorted, arbitrary order to save some run time.")
    private boolean DONT_SORT_OUTPUT = false;

    /**
     * This calculation is critically dependent on being able to skip over known polymorphic sites. Please be sure that you know what you are doing if you use this option.
     */
    @Argument(fullName="run_without_dbsnp_potentially_ruining_quality", shortName="run_without_dbsnp_potentially_ruining_quality", required=false, doc="If specified, allows the recalibrator to be used without a dbsnp rod. Very unsafe and for expert users only.")
    private boolean RUN_WITHOUT_DBSNP = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private final RecalDataManager dataManager = new RecalDataManager(); // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>(); // A list to hold the covariate objects that were requested
    private static final double DBSNP_VS_NOVEL_MISMATCH_RATE = 2.0;      // rate at which dbSNP sites (on an individual level) mismatch relative to novel sites (determined by looking at NA12878)
    private static int DBSNP_VALIDATION_CHECK_FREQUENCY = 1000000;       // how often to validate dbsnp mismatch rate (in terms of loci seen)

    public static class CountedData {
        private long countedSites = 0; // Number of loci used in the calculations, used for reporting in the output file
        private long countedBases = 0; // Number of bases used in the calculations, used for reporting in the output file
        private long skippedSites = 0; // Number of loci skipped because it was a dbSNP site, used for reporting in the output file
        private long solidInsertedReferenceBases = 0; // Number of bases where we believe SOLID has inserted the reference because the color space is inconsistent with the read base
        private long otherColorSpaceInconsistency = 0; // Number of bases where the color space is inconsistent with the read but the reference wasn't inserted.

        private long dbSNPCountsMM = 0, dbSNPCountsBases = 0;  // mismatch/base counts for dbSNP loci
        private long novelCountsMM = 0, novelCountsBases = 0;  // mismatch/base counts for non-dbSNP loci
        private int lociSinceLastDbsnpCheck = 0;               // loci since last dbsnp validation

        /**
         * Adds the values of other to this, returning this
         * @param other
         * @return this object
         */
        public CountedData add(CountedData other) {
            countedSites += other.countedSites;
            countedBases += other.countedBases;
            skippedSites += other.skippedSites;
            solidInsertedReferenceBases += other.solidInsertedReferenceBases;
            otherColorSpaceInconsistency += other.otherColorSpaceInconsistency;
            dbSNPCountsMM += other.dbSNPCountsMM;
            dbSNPCountsBases += other.dbSNPCountsBases;
            novelCountsMM += other.novelCountsMM;
            novelCountsBases += other.novelCountsBases;
            lociSinceLastDbsnpCheck += other.lociSinceLastDbsnpCheck;
            return this;
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Parse the -cov arguments and create a list of covariates to be used here
     * Based on the covariates' estimates for initial capacity allocate the data hashmap
     */
    public void initialize() {

        if( RAC.FORCE_READ_GROUP != null ) { RAC.DEFAULT_READ_GROUP = RAC.FORCE_READ_GROUP; }
        if( RAC.FORCE_PLATFORM != null ) { RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM; }

        // Get a list of all available covariates
        final List<Class<? extends Covariate>> covariateClasses = new PluginManager<Covariate>( Covariate.class ).getPlugins();
        final List<Class<? extends RequiredCovariate>> requiredClasses = new PluginManager<RequiredCovariate>( RequiredCovariate.class ).getPlugins();
        final List<Class<? extends StandardCovariate>> standardClasses = new PluginManager<StandardCovariate>( StandardCovariate.class ).getPlugins();

        // Print and exit if that's what was requested
        if ( LIST_ONLY ) {
            logger.info( "Available covariates:" );
            for( Class<?> covClass : covariateClasses ) {
                logger.info( covClass.getSimpleName() );
            }
            logger.info("");

            System.exit( 0 ); // Early exit here because user requested it
        }

        // Warn the user if no dbSNP file or other variant mask was specified
        if( knownSites.isEmpty() && !RUN_WITHOUT_DBSNP ) {
            throw new UserException.CommandLineException("This calculation is critically dependent on being able to skip over known variant sites. Please provide a VCF file containing known sites of genetic variation.");
        }

        // Initialize the requested covariates by parsing the -cov argument
        // First add the required covariates
        if( requiredClasses.size() == 2) { // readGroup and reported quality score
            requestedCovariates.add( new ReadGroupCovariate() ); // Order is important here
            requestedCovariates.add( new QualityScoreCovariate() );
        } else {
            throw new UserException.CommandLineException("There are more required covariates than expected. The instantiation list needs to be updated with the new required covariate and in the correct order.");
        }
        // Next add the standard covariates if -standard was specified by the user
        if( USE_STANDARD_COVARIATES ) {
            // We want the standard covariates to appear in a consistent order but the packageUtils method gives a random order
            // A list of Classes can't be sorted, but a list of Class names can be
            final List<String> standardClassNames = new ArrayList<String>();
            for( Class<?> covClass : standardClasses ) {
                standardClassNames.add( covClass.getName() );
            }
            Collections.sort(standardClassNames); // Sort the list of class names
            for( String className : standardClassNames ) {
                for( Class<?> covClass : standardClasses ) { // Find the class that matches this class name
                    if( covClass.getName().equals( className ) ) {
                        try {
                            final Covariate covariate = (Covariate)covClass.newInstance();
                            requestedCovariates.add( covariate );
                        } catch (Exception e) {
                            throw new DynamicClassResolutionException(covClass, e);
                        }
                    }
                }
            }
        }
        // Finally parse the -cov arguments that were provided, skipping over the ones already specified
        if( COVARIATES != null ) {
            for( String requestedCovariateString : COVARIATES ) {
                boolean foundClass = false;
                for( Class<?> covClass : covariateClasses ) {
                    if( requestedCovariateString.equalsIgnoreCase( covClass.getSimpleName() ) ) { // -cov argument matches the class name for an implementing class
                        foundClass = true;
                        if( !requiredClasses.contains( covClass ) && (!USE_STANDARD_COVARIATES || !standardClasses.contains( covClass )) ) {
                            try {
                                // Now that we've found a matching class, try to instantiate it
                                final Covariate covariate = (Covariate)covClass.newInstance();
                                requestedCovariates.add( covariate );
                            } catch (Exception e) {
                                throw new DynamicClassResolutionException(covClass, e);
                            }
                        }
                    }
                }

                if( !foundClass ) {
                    throw new UserException.CommandLineException( "The requested covariate type (" + requestedCovariateString + ") isn't a valid covariate option. Use --list to see possible covariates." );
                }
            }
        }

        logger.info( "The covariates being used here: " );
        for( Covariate cov : requestedCovariates ) {
            logger.info( "\t" + cov.getClass().getSimpleName() );
            cov.initialize( RAC ); // Initialize any covariate member variables using the shared argument collection
        }
    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * For each read at this locus get the various covariate values and increment that location in the map based on
     *   whether or not the base matches the reference at this particular location
     * @param tracker The reference metadata tracker
     * @param ref The reference context
     * @param context The alignment context
     * @return Returns 1, but this value isn't used in the reduce step
     */
    public CountedData map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        // Only use data from non-dbsnp sites
        // Assume every mismatch at a non-dbsnp site is indicative of poor quality
        CountedData counter = new CountedData();
        if( tracker.getValues(knownSites).size() == 0 ) { // If something here is in one of the knownSites tracks then skip over it, otherwise proceed
            // For each read at this locus
            for( final PileupElement p : context.getBasePileup() ) {
                final GATKSAMRecord gatkRead = (GATKSAMRecord) p.getRead();
                int offset = p.getOffset();

                if( gatkRead.containsTemporaryAttribute( SKIP_RECORD_ATTRIBUTE  ) ) {
                    continue;
                }

                if( !gatkRead.containsTemporaryAttribute( SEEN_ATTRIBUTE  ) )
                {
                    gatkRead.setTemporaryAttribute( SEEN_ATTRIBUTE, true );
                    RecalDataManager.parseSAMRecord( gatkRead, RAC );

                    // Skip over reads with no calls in the color space if the user requested it
                    if( !(RAC.SOLID_NOCALL_STRATEGY == RecalDataManager.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION) && RecalDataManager.checkNoCallColorSpace( gatkRead ) ) {
                        gatkRead.setTemporaryAttribute( SKIP_RECORD_ATTRIBUTE, true);
                        continue;
                    }

                    RecalDataManager.parseColorSpace( gatkRead );
                    gatkRead.setTemporaryAttribute( COVARS_ATTRIBUTE,
                            RecalDataManager.computeCovariates( gatkRead, requestedCovariates ));
                }


                // Skip this position if base quality is zero
                if( gatkRead.getBaseQualities()[offset] > 0 ) {

                    byte[] bases = gatkRead.getReadBases();
                    byte refBase = ref.getBase();

                    // Skip if this base is an 'N' or etc.
                    if( BaseUtils.isRegularBase( bases[offset] ) ) {

                        // SOLID bams have inserted the reference base into the read if the color space in inconsistent with the read base so skip it
                        if( !gatkRead.getReadGroup().getPlatform().toUpperCase().contains("SOLID") || RAC.SOLID_RECAL_MODE == RecalDataManager.SOLID_RECAL_MODE.DO_NOTHING ||
                            !RecalDataManager.isInconsistentColorSpace( gatkRead, offset ) ) {

                            // This base finally passed all the checks for a good base, so add it to the big data hashmap
                            updateDataFromRead( counter, gatkRead, offset, refBase );

                        } else { // calculate SOLID reference insertion rate
                            if( refBase == bases[offset] ) {
                                counter.solidInsertedReferenceBases++;
                            } else {
                                counter.otherColorSpaceInconsistency++;
                            }
                        }
                    }
                }
            }
            counter.countedSites++;
        } else { // We skipped over the dbSNP site, and we are only processing every Nth locus
            counter.skippedSites++;
            updateMismatchCounts(counter, context, ref.getBase()); // For sanity check to ensure novel mismatch rate vs dnsnp mismatch rate is reasonable
        }

        return counter;
    }

   /**
     * Update the mismatch / total_base counts for a given class of loci.
     *
     * @param counter The CountedData to be updated
     * @param context The AlignmentContext which holds the reads covered by this locus
     * @param refBase The reference base
     */
    private static void updateMismatchCounts(CountedData counter, final AlignmentContext context, final byte refBase) {
        for( PileupElement p : context.getBasePileup() ) {
            final byte readBase = p.getBase();
            final int readBaseIndex = BaseUtils.simpleBaseToBaseIndex(readBase);
            final int refBaseIndex  = BaseUtils.simpleBaseToBaseIndex(refBase);

            if( readBaseIndex != -1 && refBaseIndex != -1 ) {
                if( readBaseIndex != refBaseIndex ) {
                    counter.novelCountsMM++;
                }
                counter.novelCountsBases++;
            }
        }
    }

    /**
     * Major workhorse routine for this walker.
     * Loop through the list of requested covariates and pick out the value from the read, offset, and reference
     * Using the list of covariate values as a key, pick out the RecalDatum and increment,
     *   adding one to the number of observations and potentially one to the number of mismatches
     * Lots of things are passed as parameters to this method as a strategy for optimizing the covariate.getValue calls
     *   because pulling things out of the SAMRecord is an expensive operation.
     * @param counter Data structure which holds the counted bases
     * @param gatkRead The SAMRecord holding all the data for this read
     * @param offset The offset in the read for this locus
     * @param refBase The reference base at this locus
     */
    private void updateDataFromRead(CountedData counter, final GATKSAMRecord gatkRead, final int offset, final byte refBase) {
        final Object[][] covars = (Comparable[][]) gatkRead.getTemporaryAttribute(COVARS_ATTRIBUTE);
        final Object[] key = covars[offset];

        // Using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
        final NestedHashMap data = dataManager.data; //optimization - create local reference
        RecalDatumOptimized datum = (RecalDatumOptimized) data.get( key );
        if( datum == null ) { // key doesn't exist yet in the map so make a new bucket and add it
            // initialized with zeros, will be incremented at end of method
            datum = (RecalDatumOptimized)data.put( new RecalDatumOptimized(), true, (Object[])key );
        }

        // Need the bases to determine whether or not we have a mismatch
        final byte base = gatkRead.getReadBases()[offset];
        final long curMismatches = datum.getNumMismatches();

        // Add one to the number of observations and potentially one to the number of mismatches
        datum.incrementBaseCounts( base, refBase );
        counter.countedBases++;
        counter.novelCountsBases++;
        counter.novelCountsMM += datum.getNumMismatches() - curMismatches; // For sanity check to ensure novel mismatch rate vs dnsnp mismatch rate is reasonable
    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Initialize the reduce step by creating a PrintStream from the filename specified as an argument to the walker.
     * @return returns A PrintStream created from the -recalFile filename argument specified to the walker
     */
    public CountedData reduceInit() {
        return new CountedData();
    }

    /**
     * The Reduce method doesn't do anything for this walker.
     * @param mapped Result of the map. This value is immediately ignored.
     * @param sum The summing CountedData used to output the CSV data
     * @return returns The sum used to output the CSV data
     */
    public CountedData reduce( CountedData mapped, CountedData sum ) {
        // Do a dbSNP sanity check every so often
        return validatingDbsnpMismatchRate(sum.add(mapped));
    }

    /**
     * Validate the dbSNP reference mismatch rates.
     */
    private CountedData validatingDbsnpMismatchRate(CountedData counter) {
        if( ++counter.lociSinceLastDbsnpCheck >= DBSNP_VALIDATION_CHECK_FREQUENCY ) {
            counter.lociSinceLastDbsnpCheck = 0;

            if( counter.novelCountsBases != 0L && counter.dbSNPCountsBases != 0L ) {
                final double fractionMM_novel = (double)counter.novelCountsMM / (double)counter.novelCountsBases;
                final double fractionMM_dbsnp = (double)counter.dbSNPCountsMM / (double)counter.dbSNPCountsBases;

                if( fractionMM_dbsnp < DBSNP_VS_NOVEL_MISMATCH_RATE * fractionMM_novel ) {
                    Utils.warnUser("The variation rate at the supplied list of known variant sites seems suspiciously low. Please double-check that the correct ROD is being used. " +
                            String.format("[dbSNP variation rate = %.4f, novel variation rate = %.4f]", fractionMM_dbsnp, fractionMM_novel) );
                    DBSNP_VALIDATION_CHECK_FREQUENCY *= 2; // Don't annoyingly output the warning message every megabase of a large file
                }
            }
        }

        return counter;
    }

    public CountedData treeReduce( CountedData sum1, CountedData sum2 ) {
        return validatingDbsnpMismatchRate(sum1.add(sum2));
    }

    /**
     * Write out the full data hashmap to disk in CSV format
     * @param sum The CountedData to write out to RECAL_FILE
     */
    public void onTraversalDone( CountedData sum ) {
        logger.info( "Writing raw recalibration data..." );
        if( sum.countedBases == 0L ) {
            throw new UserException.BadInput("Could not find any usable data in the input BAM file(s).");
        }
        outputToCSV( sum, RECAL_FILE );
        logger.info( "...done!" );
    }

    /**
     * For each entry (key-value pair) in the data hashmap output the Covariate's values as well as the RecalDatum's data in CSV format
     * @param recalTableStream The PrintStream to write out to
     */
    private void outputToCSV( CountedData sum, final PrintStream recalTableStream ) {
        recalTableStream.printf("# Counted Sites    %d%n", sum.countedSites);
        recalTableStream.printf("# Counted Bases    %d%n", sum.countedBases);
        recalTableStream.printf("# Skipped Sites    %d%n", sum.skippedSites);
        recalTableStream.printf("# Fraction Skipped 1 / %.0f bp%n", (double)sum.countedSites / sum.skippedSites);

        if( sum.solidInsertedReferenceBases != 0 ) {
            recalTableStream.printf("# Fraction SOLiD inserted reference 1 / %.0f bases%n", (double) sum.countedBases / sum.solidInsertedReferenceBases);
            recalTableStream.printf("# Fraction other color space inconsistencies 1 / %.0f bases%n", (double) sum.countedBases / sum.otherColorSpaceInconsistency);
        }

        // Output header saying which covariates were used and in what order
        for( Covariate cov : requestedCovariates ) {
            recalTableStream.print( cov.getClass().getSimpleName().split("Covariate")[0] + "," );
        }
        recalTableStream.println("nObservations,nMismatches,Qempirical");

        if( DONT_SORT_OUTPUT ) {
            printMappings(recalTableStream, 0, new Object[requestedCovariates.size()], dataManager.data.data);
        } else {
            printMappingsSorted(recalTableStream, 0, new Object[requestedCovariates.size()], dataManager.data.data);
        }

        // print out an EOF marker
        recalTableStream.println(TableRecalibrationWalker.EOF_MARKER);
    }

    private void printMappingsSorted( final PrintStream recalTableStream, final int curPos, final Object[] key, final Map data) {
        final ArrayList<Comparable> keyList = new ArrayList<Comparable>();
        for( Object comp : data.keySet() ) {
            keyList.add((Comparable) comp);
        }

        Collections.sort(keyList);

        for( Comparable comp : keyList ) {
            key[curPos] = comp;
            final Object val = data.get(comp);
            if( val instanceof RecalDatumOptimized ) { // We are at the end of the nested hash maps
                // For each Covariate in the key
                for( Object compToPrint : key ) {
                    // Output the Covariate's value
                    recalTableStream.print( compToPrint + "," );
                }
                // Output the RecalDatum entry
                recalTableStream.println( ((RecalDatumOptimized)val).outputToCSV() );
            } else { // Another layer in the nested hash map
                printMappingsSorted( recalTableStream, curPos + 1, key, (Map) val );
            }
        }
    }

    private void printMappings( final PrintStream recalTableStream, final int curPos, final Object[] key, final Map data) {
        for( Object comp : data.keySet() ) {
            key[curPos] = comp;
            final Object val = data.get(comp);
            if( val instanceof RecalDatumOptimized ) { // We are at the end of the nested hash maps
                // For each Covariate in the key
                for( Object compToPrint : key ) {
                    // Output the Covariate's value
                    recalTableStream.print( compToPrint + "," );
                }
                // Output the RecalDatum entry
                recalTableStream.println( ((RecalDatumOptimized)val).outputToCSV() );
            } else { // Another layer in the nested hash map
                printMappings( recalTableStream, curPos + 1, key, (Map) val );
            }
        }
    }
}

