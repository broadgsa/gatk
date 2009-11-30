package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.ArgumentCollection;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.util.*;

import net.sf.samtools.SAMRecord;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 *
 * This walker is designed to work as the first pass in a two-pass processing step.
 * It does a by-locus traversal operating only at sites that are not in dbSNP.
 * We assume that all reference mismatches we see are therefore errors and indicitive of poor base quality.
 * This walker generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and dinucleotide)
 * Since there is a large amount of data one can then calculate an empirical probability of error
 *   given the particular covariates seen at this site, where p(error) = num mismatches / num observations
 * The output file is a CSV list of (the several covariate values, num observations, num mismatches, empirical quality score)
 * The first lines of the output file give the name of the covariate classes that were used for this calculation.
 *
 * Note: ReadGroupCovariate and QualityScoreCovariate are required covariates and must be at the start of the list.
 * Note: This walker is designed to be used in conjunction with TableRecalibrationWalker.
 */

@By( DataSource.READS ) // Only look at covered loci, not every loci of the reference file
@WalkerName( "CountCovariates" )
@ReadFilters( {ZeroMappingQualityReadFilter.class} ) // Filter out all reads with zero mapping quality
@Requires( {DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES} ) // This walker requires both -I input.bam and -R reference.fasta
public class CovariateCounterWalker extends LocusWalker<Integer, PrintStream> {

    /////////////////////////////
    // Shared Arguments
    /////////////////////////////
    @ArgumentCollection private RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="list", shortName="ls", doc="List the available covariates and exit", required=false)
    private Boolean LIST_ONLY = false;
    @Argument(fullName="covariate", shortName="cov", doc="Covariates to be used in the recalibration. Each covariate is given as a separate cov parameter. ReadGroup and ReportedQuality are already added for you.", required=false)
    private String[] COVARIATES = null;
    @Argument(fullName="process_nth_locus", shortName="pN", required=false, doc="Only process every Nth covered locus we see.")
    private int PROCESS_EVERY_NTH_LOCUS = 1;

    /////////////////////////////
    // Debugging-only Arguments
    /////////////////////////////
    @Argument(fullName="no_print_header", shortName="noHeader", required=false, doc="Don't print the usual header on the table recalibration file. FOR DEBUGGING PURPOSES ONLY.")
    private boolean NO_PRINT_HEADER = false;
    @Argument(fullName="sorted_output", shortName="sorted", required=false, doc="The outputted table recalibration file will be in sorted order at the cost of added overhead. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.")
    private boolean SORTED_OUTPUT = false;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private RecalDataManager dataManager; // Holds the data HashMap, mostly used by TableRecalibrationWalker to create collapsed data hashmaps
    private ArrayList<Covariate> requestedCovariates; // A list to hold the covariate objects that were requested
    private IdentityHashMap<SAMRecord, ReadHashDatum> readDatumHashMap; // A hash map that hashes the read object itself into properties commonly pulled out of the read. Done for optimization purposes.
    private int sizeOfReadDatumHashMap = 0;
    private long countedSites = 0; // Number of loci used in the calculations, used for reporting in the output file
    private long countedBases = 0; // Number of bases used in the calculations, used for reporting in the output file
    private long skippedSites = 0; // Number of loci skipped because it was a dbSNP site, used for reporting in the output file
    private int numUnprocessed = 0; // Number of consecutive loci skipped because we are only processing every Nth site
    private static final String versionString = "v2.0.8"; // Major version, minor version, and build number
    private Pair<Long, Long> dbSNP_counts = new Pair<Long, Long>(0L, 0L);  // mismatch/base counts for dbSNP loci
    private Pair<Long, Long> novel_counts = new Pair<Long, Long>(0L, 0L);  // mismatch/base counts for non-dbSNP loci
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

        logger.info( "CovariateCounterWalker version: " + versionString );
        if( RAC.FORCE_READ_GROUP != null ) { RAC.DEFAULT_READ_GROUP = RAC.FORCE_READ_GROUP; }
        if( RAC.FORCE_PLATFORM != null ) { RAC.DEFAULT_PLATFORM = RAC.FORCE_PLATFORM; }

        // Get a list of all available covariates
        final List<Class<? extends Covariate>> classes = PackageUtils.getClassesImplementingInterface( Covariate.class );
   
        // Print and exit if that's what was requested
        if ( LIST_ONLY ) {
            out.println( "Available covariates:" );
            for( Class<?> covClass : classes ) {
                out.println( covClass.getSimpleName() );
            }
            out.println();
            
            System.exit( 0 ); // Early exit here because user requested it
        }

        // Warn the user if no dbSNP file was specified
        boolean foundDBSNP = false;
        for( ReferenceOrderedDataSource rod : this.getToolkit().getRodDataSources() ) {
            if( rod.getName().equalsIgnoreCase( "dbsnp" ) ) {
                foundDBSNP = true;
            }
        }
        if( !foundDBSNP ) {
            Utils.warnUser("This calculation is critically dependent on being able to skip over known variant sites. Are you sure you want to be running without a dbSNP rod specified?");
        }


        // Initialize the requested covariates by parsing the -cov argument
        // BUGBUG: This is a mess because there are a lot of cases (validate, all, none, and supplied covList). Clean up needed.
        requestedCovariates = new ArrayList<Covariate>();
        int estimatedCapacity = 1; // Capacity is multiplicitive so this starts at one
        if( RAC.VALIDATE_OLD_RECALIBRATOR ) {
            requestedCovariates.add( new ReadGroupCovariate() );
            requestedCovariates.add( new CycleCovariate() ); // Unfortunately order is different here in order to match the old recalibrator exactly
            requestedCovariates.add( new QualityScoreCovariate() );
            requestedCovariates.add( new DinucCovariate() );
            estimatedCapacity = 60 * 100 * 40 * 16;
        } else if( COVARIATES != null ) {
            if(COVARIATES[0].equalsIgnoreCase( "ALL" )) { // The user wants ALL covariates to be used
                requestedCovariates.add( new ReadGroupCovariate() ); // First add the required covariates then add the rest by looping over all implementing classes that were found
                requestedCovariates.add( new QualityScoreCovariate() );
                for( Class<?> covClass : classes ) {
                    try {
                        Covariate covariate = (Covariate)covClass.newInstance();
                        
                        estimatedCapacity *= covariate.estimatedNumberOfBins();
                        if( !( covariate instanceof ReadGroupCovariate || covariate instanceof QualityScoreCovariate ) ) { // These were already added so don't add them again
                            requestedCovariates.add( covariate );
                        }
                    } catch ( InstantiationException e ) {
                        throw new StingException( String.format("Can not instantiate covariate class '%s': must be concrete class.", covClass.getSimpleName()) );
                    } catch ( IllegalAccessException e ) {
                        throw new StingException( String.format("Can not instantiate covariate class '%s': must have no-arg constructor.", covClass.getSimpleName()) );
                    }
               }
            } else { // The user has specified a list of several covariates
                int covNumber = 1;
                for( String requestedCovariateString : COVARIATES ) {
                    boolean foundClass = false;
                    for( Class<?> covClass : classes ) {
                        if( requestedCovariateString.equalsIgnoreCase( covClass.getSimpleName() ) ) { // -cov argument matches the class name for an implementing class
                            foundClass = true;
                            // Read Group Covariate and Quality Score Covariate are required covariates for the recalibration calculation and must begin the list
                            if( (covNumber == 1 && !requestedCovariateString.equalsIgnoreCase( "ReadGroupCovariate" )) ||
                                (covNumber == 2 && !requestedCovariateString.equalsIgnoreCase( "QualityScoreCovariate" )) ) {
                                throw new StingException("ReadGroupCovariate and QualityScoreCovariate are required covariates for the recalibration calculation and must begin the list" );
                            }
                            covNumber++;
                            try {
                                // Now that we've found a matching class, try to instantiate it
                                Covariate covariate = (Covariate)covClass.newInstance();
                                estimatedCapacity *= covariate.estimatedNumberOfBins();
                                requestedCovariates.add( covariate );
                            } catch ( InstantiationException e ) {
                                throw new StingException( String.format("Can not instantiate covariate class '%s': must be concrete class.", covClass.getSimpleName()) );
                            } catch ( IllegalAccessException e ) {
                                throw new StingException( String.format("Can not instantiate covariate class '%s': must have no-arg constructor.", covClass.getSimpleName()) );
                            }
                        }
                    }

                    if( !foundClass ) {
                        throw new StingException( "The requested covariate type (" + requestedCovariateString + ") isn't a valid covariate option. Use --list to see possible covariates." );
                    }
                }
            }
        } else { // No covariates were specified by the user so add the default, required ones
            Utils.warnUser( "Using default set of covariates because none were specified. Using ReadGroupCovariate and QualityScoreCovariate only." );
            requestedCovariates.add( new ReadGroupCovariate() );
            requestedCovariates.add( new QualityScoreCovariate() );
            estimatedCapacity = 60 * 40;
        }

        logger.info( "The covariates being used here: " );
        for( Covariate cov : requestedCovariates ) {
            logger.info( "\t" + cov.getClass().getSimpleName() );
            cov.initialize( RAC ); // Initialize any covariate member variables using the shared argument collection
        }

        // Don't want to crash with out of heap space exception
        if( estimatedCapacity > 300 * 40 * 200 || estimatedCapacity < 0 ) { // Could be negative if overflowed
            estimatedCapacity = 300 * 40 * 200;
        }
        dataManager = new RecalDataManager( estimatedCapacity );
        readDatumHashMap = new IdentityHashMap<SAMRecord, ReadHashDatum>();
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

        // Pull out anything passed by -B name,type,file that has the name "dbsnp"
        final RODRecordList<ReferenceOrderedDatum> dbsnpRODs = tracker.getTrackData( "dbsnp", null );
        boolean isSNP = false;
        if (dbsnpRODs != null) {
            for( ReferenceOrderedDatum rod : dbsnpRODs ) {
                if( ((Variation)rod).isSNP() ) {
                    isSNP = true; // At least one of the rods says this is a snp site
                    break;
                }
            }
        }

        // Only use data from non-dbsnp sites
        // Assume every mismatch at a non-dbsnp site is indicitive of poor quality
        if( !isSNP && ( ++numUnprocessed >= PROCESS_EVERY_NTH_LOCUS ) ) {
            numUnprocessed = 0; // Reset the counter because we are processing this very locus

            // Long list of variables that are preallocated for use in for loop below
            final List<SAMRecord> reads = context.getReads();
            final List<Integer> offsets = context.getOffsets();
            SAMRecord read;
            int offset;
            byte refBase;
            byte prevBase;
            byte[] colorSpaceQuals;
            ReadHashDatum readDatum;

            final int numReads = reads.size();
            // For each read at this locus
            for( int iii = 0; iii < numReads; iii++ ) {
                read = reads.get(iii);
                offset = offsets.get(iii); // offset is zero based so quals[offset] and bases[offset] is correct

                // Try to pull the read out of the read IdentityHashMap
                readDatum = readDatumHashMap.get( read );
                if( readDatum == null ) {

                    // If the HashMap of read objects has grown too large then throw out the (mostly stale) reads
                    if( sizeOfReadDatumHashMap > 100000 ) { //BUGBUG: Can I make this number larger?
                        readDatumHashMap.clear();
                        sizeOfReadDatumHashMap = 0;
                    }

                    // This read isn't in the hashMap yet so fill out the datum and add it to the map so that we never have to do the work again
                    readDatum = ReadHashDatum.parseSAMRecord( read, RAC );
                    readDatumHashMap.put( read, readDatum );
                    sizeOfReadDatumHashMap++;
                }

                // Skip first and last base because there is no dinuc
                // BUGBUG: Technically we only have to skip the first base on forward reads and the last base on negative strand reads. Change after passing old integration tests.
                if( offset > 0 ) {
                    if( offset < readDatum.length - 1 ) {
                        // Skip if base quality is zero
                        if( readDatum.quals[offset] > 0 ) {

                            refBase = (byte)ref.getBase();
                            prevBase = readDatum.bases[offset - 1];

                            // DinucCovariate is responsible for getting the complement bases if needed
                            if( readDatum.isNegStrand ) {
                                prevBase = readDatum.bases[offset + 1];
                            }

                            // Skip if this base or the previous one was an 'N' or etc.
                            if( BaseUtils.isRegularBase( (char)prevBase ) && BaseUtils.isRegularBase( (char)(readDatum.bases[offset]) ) ) {

                                // SOLID bams insert the reference base into the read if the color space quality is zero, so skip over them
                                colorSpaceQuals = null;
                                if( readDatum.platform.equalsIgnoreCase("SOLID") ) {
                                    colorSpaceQuals = QualityUtils.fastqToPhred((String)read.getAttribute(RecalDataManager.COLOR_SPACE_QUAL_ATTRIBUTE_TAG));
                                }
                                if( colorSpaceQuals == null || colorSpaceQuals[offset] > 0 ) //BUGBUG: This isn't exactly correct yet
                                {
                                    // This base finally passed all the checks for a good base, so add it to the big data hashmap
                                    updateDataFromRead( readDatum, offset, refBase );
                                }
                            } else {
                                if( RAC.VALIDATE_OLD_RECALIBRATOR ) {
                                    countedBases++; // Replicating a small bug in the old recalibrator
                                }
                            }
                        }
                    } else { // At the last base in the read so we can remove the read from our IdentityHashMap since we will never see it again
                        readDatumHashMap.remove( read );
                        sizeOfReadDatumHashMap--;
                    }
                }
            }
            countedSites++;
            updateMismatchCounts(novel_counts, context, ref.getBase()); // For sanity check to ensure novel mismatch rate vs dnsnp mismatch rate is reasonable

        } else { // We skipped over the dbSNP site, and we are only processing every Nth locus
            skippedSites++;
            if( isSNP) {
                updateMismatchCounts(dbSNP_counts, context, ref.getBase());// For sanity check to ensure novel mismatch rate vs dnsnp mismatch rate is reasonable
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
    private static void updateMismatchCounts(Pair<Long, Long> counts, AlignmentContext context, char ref) {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        for(int iii = 0; iii < reads.size(); iii++ ) {
            char readChar = (char)(reads.get(iii).getReadBases()[offsets.get(iii)]);
            int readCharBaseIndex = BaseUtils.simpleBaseToBaseIndex(readChar);
            int refCharBaseIndex  = BaseUtils.simpleBaseToBaseIndex(ref);

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
        if( novel_counts.second == 0 || dbSNP_counts.second == 0 ) {
            return;
        }

        double fractionMM_novel = (double)novel_counts.first / (double)novel_counts.second;
        double fractionMM_dbsnp = (double)dbSNP_counts.first / (double)dbSNP_counts.second;
        
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
     * @param readDatum The ReadHashDatum holding all the important properties of this read
     * @param offset The offset in the read for this locus
     * @param refBase The reference base at this locus
     */
    private void updateDataFromRead(final ReadHashDatum readDatum, final int offset, final byte refBase) {

        List<Comparable> key = new ArrayList<Comparable>();
        
        // Loop through the list of requested covariates and pick out the value from the read, offset, and reference
        for( Covariate covariate : requestedCovariates ) {
        	key.add( covariate.getValue( readDatum, offset ) );
        }

    	// Using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
        RecalDatum datum = dataManager.data.get( key );
        if( datum == null ) { // key doesn't exist yet in the map so make a new bucket and add it
            datum = new RecalDatum(); // initialized with zeros, will be incremented at end of method
            if( RAC.VALIDATE_OLD_RECALIBRATOR || SORTED_OUTPUT ) {
                dataManager.data.sortedPut( key, datum );
            } else {
                dataManager.data.put( key, datum );
            }
        }
        
        // Need the bases to determine whether or not we have a mismatch
        byte base = readDatum.bases[offset];
        
        // Add one to the number of observations and potentially one to the number of mismatches
        datum.increment( (char)base, (char)refBase ); // Dangerous: If you don't cast to char than the bytes default to the (long, long) version of the increment method which is really bad
        countedBases++;
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

        //BUGBUG: This method is a mess. It will be cleaned up when I get rid of the validation and no_header debug options.
        if( RAC.VALIDATE_OLD_RECALIBRATOR ) {
            // Output the old header as well as output the data in sorted order
            recalTableStream.printf("# collapsed_pos    false%n");
            recalTableStream.printf("# collapsed_dinuc  false%n");
            recalTableStream.printf("# counted_sites    %d%n", countedSites);
            recalTableStream.printf("# counted_bases    %d%n", countedBases);
            recalTableStream.printf("# skipped_sites    %d%n", skippedSites);
            recalTableStream.printf("# fraction_skipped 1 / %.0f bp%n", (double)countedSites / skippedSites);
            recalTableStream.printf("rg,pos,Qrep,dn,nBases,nMismatches,Qemp%n");
            for( Pair<List<? extends Comparable>,RecalDatum> entry : dataManager.data.entrySetSorted4() ) {
                for( Comparable comp : entry.first ) {
                    recalTableStream.print( comp + "," );
                }
                recalTableStream.println( entry.second.outputToCSV() );
            }
        } else {
            if( !NO_PRINT_HEADER ) {
                recalTableStream.printf("# Counted Sites    %d%n", countedSites);
                recalTableStream.printf("# Counted Bases    %d%n", countedBases);
                recalTableStream.printf("# Skipped Sites    %d%n", skippedSites);
                if( PROCESS_EVERY_NTH_LOCUS == 1 ) {
                    recalTableStream.printf("# Fraction Skipped 1 / %.0f bp%n", (double)countedSites / skippedSites);
                } else {
                    recalTableStream.printf("# Percent Skipped  %.4f%n", 100.0 * (double)skippedSites / ((double)countedSites+skippedSites));
                }
                for( Covariate cov : requestedCovariates ) {
                    // The "@!" is a code for TableRecalibrationWalker to recognize this line as a Covariate class name
                    recalTableStream.print( cov.getClass().getSimpleName().split("Covariate")[0] + "," );
                }
                recalTableStream.println("nObservations,nMismatches,Qempirical");
            }

            if( SORTED_OUTPUT && requestedCovariates.size() == 4 )
            {
                for( Pair<List<? extends Comparable>,RecalDatum> entry : dataManager.data.entrySetSorted4() ) { //BUGBUG: entrySetSorted4 isn't correct here
                    for( Comparable comp : entry.first ) {
                        recalTableStream.print( comp + "," );
                    }
                    recalTableStream.println( entry.second.outputToCSV() );
                }
            } else {
                // For each entry in the data hashmap
                for( Map.Entry<List<? extends Comparable>, RecalDatum> entry : dataManager.data.entrySet() ) {
                    // For each Covariate in the key
                    for( Comparable comp : entry.getKey() ) {
                        // Output the Covariate's value
                        recalTableStream.print( comp + "," );
                    }
                    // Output the RecalDatum entry
                    recalTableStream.println( entry.getValue().outputToCSV() );
                }
            }
        }
    }
}

