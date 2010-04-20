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

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.classloader.PackageUtils;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * This walker is designed to work as the first pass in a two-pass processing step.
 * It does a by-locus traversal operating only at sites that are not in dbSNP.
 * We assume that all reference mismatches we see are therefore errors and indicative of poor base quality.
 * This walker generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and dinucleotide)
 * Since there is a large amount of data one can then calculate an empirical probability of error
 *   given the particular covariates seen at this site, where p(error) = num mismatches / num observations
 * The output file is a CSV list of (the several covariate values, num observations, num mismatches, empirical quality score)
 * The first non-comment line of the output file gives the name of the covariates that were used for this calculation.
 *
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and will be added for the user regardless of whether or not they were specified
 * Note: This walker is designed to be used in conjunction with TableRecalibrationWalker.
 *
 * @author rpoplin
 * @since Nov 3, 2009
 * @help.summary First pass of the recalibration. Generates recalibration table based on various user-specified covariates (such as reported quality score, cycle, and dinucleotide).
 */

@By( DataSource.READS ) // Only look at covered loci, not every loci of the reference file
@WalkerName( "CountCovariates" )
@ReadFilters( {ZeroMappingQualityReadFilter.class} ) // Filter out all reads with zero mapping quality
@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta
public class CovariateCounterWalker extends LocusWalker<Integer, PrintStream> {

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
    @Argument(fullName="list", shortName="ls", doc="List the available covariates and exit", required=false)
    private boolean LIST_ONLY = false;
    @Argument(fullName="covariate", shortName="cov", doc="Covariates to be used in the recalibration. Each covariate is given as a separate cov parameter. ReadGroup and ReportedQuality are required covariates and are already added for you.", required=false)
    private String[] COVARIATES = null;
    @Argument(fullName="standard_covs", shortName="standard", doc="Use the standard set of covariates in addition to the ones listed using the -cov argument", required=false)
    private boolean USE_STANDARD_COVARIATES = false;
    @Argument(fullName="process_nth_locus", shortName="pN", required=false, doc="Only process every Nth covered locus we see.")
    private int PROCESS_EVERY_NTH_LOCUS = 1;

    /////////////////////////////
    // Debugging-only Arguments
    /////////////////////////////
    @Argument(fullName="dont_sort_output", shortName="unsorted", required=false, doc="If specified, the output table recalibration csv file will be in an unsorted, arbitrary order to save some run time.")
    private boolean DONT_SORT_OUTPUT = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private final RecalDataManager dataManager = new RecalDataManager(); // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private final ArrayList<Covariate> requestedCovariates = new ArrayList<Covariate>(); // A list to hold the covariate objects that were requested
    private long countedSites = 0; // Number of loci used in the calculations, used for reporting in the output file
    private long countedBases = 0; // Number of bases used in the calculations, used for reporting in the output file
    private long skippedSites = 0; // Number of loci skipped because it was a dbSNP site, used for reporting in the output file
    private long solidInsertedReferenceBases = 0; // Number of bases where we believe SOLID has inserted the reference because the color space is inconsistent with the read base
    private long otherColorSpaceInconsistency = 0; // Number of bases where the color space is inconsistent with the read but the reference wasn't inserted.
    private int numUnprocessed = 0; // Number of consecutive loci skipped because we are only processing every Nth site
    private final Pair<Long, Long> dbSNP_counts = new Pair<Long, Long>(0L, 0L);  // mismatch/base counts for dbSNP loci
    private final Pair<Long, Long> novel_counts = new Pair<Long, Long>(0L, 0L);  // mismatch/base counts for non-dbSNP loci
    private static final double DBSNP_VS_NOVEL_MISMATCH_RATE = 2.0;        // rate at which dbSNP sites (on an individual level) mismatch relative to novel sites (determined by looking at NA12878)
    private int DBSNP_VALIDATION_CHECK_FREQUENCY = 1000000;                // how often to validate dbsnp mismatch rate (in terms of loci seen)
    private int lociSinceLastDbsnpCheck = 0;                               // loci since last dbsnp validation

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
        DBSNP_VALIDATION_CHECK_FREQUENCY *= PROCESS_EVERY_NTH_LOCUS;
        if( !RAC.checkSolidRecalMode() ) {
            throw new StingException( "Unrecognized --solid_recal_mode argument. Implemented options: DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS");
        }

        // Get a list of all available covariates
        final List<Class<? extends Covariate>> covariateClasses = PackageUtils.getClassesImplementingInterface( Covariate.class );
        final List<Class<? extends RequiredCovariate>> requiredClasses = PackageUtils.getClassesImplementingInterface( RequiredCovariate.class );
        final List<Class<? extends StandardCovariate>> standardClasses = PackageUtils.getClassesImplementingInterface( StandardCovariate.class );

        // Print and exit if that's what was requested
        if ( LIST_ONLY ) {
            out.println( "Available covariates:" );
            for( Class<?> covClass : covariateClasses ) {
                out.println( covClass.getSimpleName() );
            }
            out.println();

            System.exit( 0 ); // Early exit here because user requested it
        }

        // Warn the user if no dbSNP file was specified
        boolean foundDBSNP = false;
        for( ReferenceOrderedDataSource rod : this.getToolkit().getRodDataSources() ) {
            if( rod != null ) {
                foundDBSNP = true;
                break;
            }
        }
        if( !foundDBSNP ) {
            Utils.warnUser("This calculation is critically dependent on being able to skip over known variant sites. Are you sure you want to be running without a dbSNP rod specified?");
        }

        // Initialize the requested covariates by parsing the -cov argument
        // First add the required covariates
        if( requiredClasses.size() == 2) { // readGroup and reported quality score
            requestedCovariates.add( new ReadGroupCovariate() ); // Order is important here
            requestedCovariates.add( new QualityScoreCovariate() );
        } else {
            throw new StingException("There are more required covariates than expected. The instantiation list needs to be updated with the new required covariate and in the correct order.");
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
                        } catch ( InstantiationException e ) {
                            throw new StingException( String.format("Can not instantiate covariate class '%s': must be concrete class.", covClass.getSimpleName()) );
                        } catch ( IllegalAccessException e ) {
                            throw new StingException( String.format("Can not instantiate covariate class '%s': must have no-arg constructor.", covClass.getSimpleName()) );
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
                            } catch ( InstantiationException e ) {
                                throw new StingException( String.format("Can not instantiate covariate class '%s': must be concrete class.", covClass.getSimpleName()) );
                            } catch ( IllegalAccessException e ) {
                                throw new StingException( String.format("Can not instantiate covariate class '%s': must have no-arg constructor.", covClass.getSimpleName()) );
                            }
                        }
                    }
                }

                if( !foundClass ) {
                    throw new StingException( "The requested covariate type (" + requestedCovariateString + ") isn't a valid covariate option. Use --list to see possible covariates." );
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
    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        // Pull out data for this locus for all the input RODs and check if this is a known variant site in any of them
        boolean isSNP = false;
        for( final VariantContext vc : tracker.getAllVariantContexts(ref, null, context.getLocation(), false, false) ) {
            if( vc != null && vc.isSNP() ) {
                isSNP = true;
                break;
            }
        }
        
        // Only use data from non-dbsnp sites
        // Assume every mismatch at a non-dbsnp site is indicative of poor quality
        if( !isSNP && ( ++numUnprocessed >= PROCESS_EVERY_NTH_LOCUS ) ) {
            numUnprocessed = 0; // Reset the counter because we are processing this very locus

            GATKSAMRecord gatkRead;
            int offset;
            byte refBase;
            byte[] bases;

            // For each read at this locus
            for( PileupElement p : context.getBasePileup() ) {
                gatkRead = (GATKSAMRecord) p.getRead();
                offset = p.getOffset();

                if( gatkRead.containsTemporaryAttribute( SKIP_RECORD_ATTRIBUTE  ) ) {
                    continue;
                }

                if( !gatkRead.containsTemporaryAttribute( SEEN_ATTRIBUTE  ) )
                {
                    gatkRead.setTemporaryAttribute( SEEN_ATTRIBUTE, true );
                    RecalDataManager.parseSAMRecord( gatkRead, RAC );

                    // Skip over reads with no calls in the color space if the user requested it
                    if( RAC.IGNORE_NOCALL_COLORSPACE && RecalDataManager.checkNoCallColorSpace( gatkRead ) ) {
                        gatkRead.setTemporaryAttribute( SKIP_RECORD_ATTRIBUTE, true);
                        continue;
                    }

                    RecalDataManager.parseColorSpace( gatkRead );
                    gatkRead.setTemporaryAttribute( COVARS_ATTRIBUTE,
                            RecalDataManager.computeCovariates( gatkRead, requestedCovariates ));
                }


                // Skip this position if base quality is zero
                if( gatkRead.getBaseQualities()[offset] > 0 ) {

                    bases = gatkRead.getReadBases();
                    refBase = (byte)ref.getBase();

                    // Skip if this base is an 'N' or etc.
                    if( BaseUtils.isRegularBase( (char)(bases[offset]) ) ) {

                        // SOLID bams have inserted the reference base into the read if the color space in inconsistent with the read base so skip it
                        if( !gatkRead.getReadGroup().getPlatform().toUpperCase().contains("SOLID") || RAC.SOLID_RECAL_MODE.equalsIgnoreCase("DO_NOTHING") || !RecalDataManager.isInconsistentColorSpace( gatkRead, offset ) ) {

                            // This base finally passed all the checks for a good base, so add it to the big data hashmap
                            updateDataFromRead( gatkRead, offset, refBase );

                        } else { // calculate SOLID reference insertion rate
                            if( refBase == (char)bases[offset] ) {
                                solidInsertedReferenceBases++;
                            } else {
                                otherColorSpaceInconsistency++;
                            }
                        }
                    }
                }
            }
            countedSites++;
        } else { // We skipped over the dbSNP site, and we are only processing every Nth locus
            skippedSites++;
            if( isSNP ) {
                updateMismatchCounts(dbSNP_counts, context, ref.getBase()); // For sanity check to ensure novel mismatch rate vs dnsnp mismatch rate is reasonable
            }
        }

        // Do a dbSNP sanity check every so often
        if( ++lociSinceLastDbsnpCheck == DBSNP_VALIDATION_CHECK_FREQUENCY ) {
            lociSinceLastDbsnpCheck = 0;
            validateDbsnpMismatchRate();
        }

        return 1; // This value isn't actually used anywhere
    }



   /**
     * Update the mismatch / total_base counts for a given class of loci.
     *
     * @param counts The counts to be updated
     * @param context The AlignmentContext which holds the reads covered by this locus
     * @param ref The reference base
     */
    private static void updateMismatchCounts(final Pair<Long, Long> counts, final AlignmentContext context, final char ref) {
        for( PileupElement p : context.getBasePileup() ) {
            final char readChar = (char)(p.getBase());
            final int readCharBaseIndex = BaseUtils.simpleBaseToBaseIndex(readChar);
            final int refCharBaseIndex  = BaseUtils.simpleBaseToBaseIndex(ref);

            if( readCharBaseIndex != -1 && refCharBaseIndex != -1 ) {
                if( readCharBaseIndex != refCharBaseIndex ) {
                    counts.first++;
                }
                counts.second++;
            }
        }
    }

   /**
     * Validate the dbSNP reference mismatch rates.
     */
    private void validateDbsnpMismatchRate() {
        if( novel_counts.second == 0L || dbSNP_counts.second == 0L ) {
            return;
        }

        final double fractionMM_novel = (double)novel_counts.first / (double)novel_counts.second;
        final double fractionMM_dbsnp = (double)dbSNP_counts.first / (double)dbSNP_counts.second;

        if( fractionMM_dbsnp < DBSNP_VS_NOVEL_MISMATCH_RATE * fractionMM_novel ) {
            Utils.warnUser("The variation rate at the supplied list of known variant sites seems suspiciously low. Please double-check that the correct ROD is being used. " +
                            String.format("[dbSNP variation rate = %.4f, novel variation rate = %.4f]", fractionMM_dbsnp, fractionMM_novel) );
            DBSNP_VALIDATION_CHECK_FREQUENCY *= 2; // Don't annoyingly output the warning message every megabase of a large file
        }
    }

    /**
     * Major workhorse routine for this walker.
     * Loop through the list of requested covariates and pick out the value from the read, offset, and reference
     * Using the list of covariate values as a key, pick out the RecalDatum and increment,
     *   adding one to the number of observations and potentially one to the number of mismatches
     * Lots of things are passed as parameters to this method as a strategy for optimizing the covariate.getValue calls
     *   because pulling things out of the SAMRecord is an expensive operation.
     * @param gatkRead The SAMRecord holding all the data for this read
     * @param offset The offset in the read for this locus
     * @param refBase The reference base at this locus
     */
    private void updateDataFromRead(final GATKSAMRecord gatkRead, final int offset, final byte refBase) {


        final Object[][] covars = (Comparable[][]) gatkRead.getTemporaryAttribute(COVARS_ATTRIBUTE);
        final Object[] key = covars[offset];

        // Using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
        final NestedHashMap data = dataManager.data; //optimization - create local reference
        RecalDatumOptimized datum = (RecalDatumOptimized) data.get( key );
        if( datum == null ) { // key doesn't exist yet in the map so make a new bucket and add it
            datum = new RecalDatumOptimized(); // initialized with zeros, will be incremented at end of method
            data.put( datum, (Object[])key );
        }

        // Need the bases to determine whether or not we have a mismatch
        final byte base = gatkRead.getReadBases()[offset];
        final long curMismatches = datum.getNumMismatches();

        // Add one to the number of observations and potentially one to the number of mismatches
        datum.increment( (char)base, (char)refBase ); // Dangerous: If you don't cast to char then the bytes default to the (long, long) version of the increment method which is really bad
        countedBases++;
        novel_counts.second++;
        novel_counts.first += datum.getNumMismatches() - curMismatches; // For sanity check to ensure novel mismatch rate vs dnsnp mismatch rate is reasonable
    }


    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Initialize the reudce step by creating a PrintStream from the filename specified as an argument to the walker.
     * @return returns A PrintStream created from the -recalFile filename argument specified to the walker
     */
    public PrintStream reduceInit() {
        try {
            return new PrintStream( RAC.RECAL_FILE );
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException( "Couldn't open output file: ", e );
        }
    }

    /**
     * The Reduce method doesn't do anything for this walker.
     * @param value Result of the map. This value is immediately ignored.
     * @param recalTableStream The PrintStream used to output the CSV data
     * @return returns The PrintStream used to output the CSV data
     */
    public PrintStream reduce( Integer value, PrintStream recalTableStream ) {
        return recalTableStream; // Nothing to do here
    }

    /**
     * Write out the full data hashmap to disk in CSV format
     * @param recalTableStream The PrintStream to write out to
     */
    public void onTraversalDone( PrintStream recalTableStream ) {
        logger.info( "Writing raw recalibration data..." );
        outputToCSV( recalTableStream );
        logger.info( "...done!" );

        recalTableStream.close();
    }

    /**
     * For each entry (key-value pair) in the data hashmap output the Covariate's values as well as the RecalDatum's data in CSV format
     * @param recalTableStream The PrintStream to write out to
     */
    private void outputToCSV( final PrintStream recalTableStream ) {

        recalTableStream.printf("# Counted Sites    %d%n", countedSites);
        recalTableStream.printf("# Counted Bases    %d%n", countedBases);
        recalTableStream.printf("# Skipped Sites    %d%n", skippedSites);
        if( PROCESS_EVERY_NTH_LOCUS == 1 ) {
            recalTableStream.printf("# Fraction Skipped 1 / %.0f bp%n", (double)countedSites / skippedSites);
        } else {
            recalTableStream.printf("# Percent Skipped  %.4f%n", 100.0 * (double)skippedSites / ((double)countedSites+skippedSites));
        }
        if( solidInsertedReferenceBases != 0 ) {
            recalTableStream.printf("# Fraction SOLiD inserted reference 1 / %.0f bases%n", (double) countedBases / solidInsertedReferenceBases);
            recalTableStream.printf("# Fraction other color space inconsistencies 1 / %.0f bases%n", (double) countedBases / otherColorSpaceInconsistency);
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

