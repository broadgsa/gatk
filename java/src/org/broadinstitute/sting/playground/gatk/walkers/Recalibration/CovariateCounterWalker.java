package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.PackageUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.util.*;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 */

@WalkerName("CountCovariatesRefactored")
public class CovariateCounterWalker extends LocusWalker<Integer, PrintStream> {

    @Argument(fullName="list", shortName="ls", doc="List the available covariates and exit", required=false)
    protected Boolean LIST_ONLY = false;
    @Argument(fullName="covariate", shortName="cov", doc="Covariates to be used in the recalibration. Each covariate is given as a separate cov parameter. ReadGroup and ReportedQuality are already added for you.", required=false)
    protected String[] COVARIATES = null;
    @Argument(fullName="min_mapping_quality", shortName="minmap", required=false, doc="Only use reads with at least this mapping quality score")
    public int MIN_MAPPING_QUALITY = 1;
    @Argument(fullName = "use_original_quals", shortName="OQ", doc="If provided, we will use use the quals from the original qualities OQ attribute field instead of the quals in the regular QUALS field", required=false)
    public boolean USE_ORIGINAL_QUALS = false;
    @Argument(fullName="recal_file", shortName="rf", required=false, doc="Filename for the outputted covariates table recalibration file")
    public String RECAL_FILE = "output.recal_data.csv";

    protected static RecalDataManager dataManager;
    protected static ArrayList<Covariate> requestedCovariates;

    public void initialize() {

        dataManager = new RecalDataManager();

        // Get a list of all available covariates
        List<Class<? extends Covariate>> classes = PackageUtils.getClassesImplementingInterface(Covariate.class);

        // Print and exit if that's what was requested
        if ( LIST_ONLY ) {
            out.println( "Available covariates:" );
            for( Class<?> covClass : classes ) {
                out.println( covClass.getSimpleName() );
            }
            out.println();
            
            System.exit( 0 ); // early exit here because user requested it
        }

        // Initialize the requested covariates
        requestedCovariates = new ArrayList<Covariate>();
        requestedCovariates.add( new ReadGroupCovariate() );    // Read Group Covariate is a required covariate for the recalibration calculation
        requestedCovariates.add( new QualityScoreCovariate() ); // Quality Score Covariate is a required covariate for the recalibration calculation
        if( COVARIATES != null ) {
            for( String requestedCovariateString : COVARIATES ) {

                boolean foundClass = false;
                for( Class<?> covClass : classes ) {

                    if( requestedCovariateString.equalsIgnoreCase( covClass.getSimpleName() ) ) {
                        foundClass = true;
                        try {
                            Covariate covariate = (Covariate)covClass.newInstance();
                            requestedCovariates.add( covariate );
                            if (covariate instanceof ReadGroupCovariate || covariate instanceof QualityScoreCovariate) {
                                throw new StingException( "ReadGroupCovariate and QualityScoreCovariate are required covariates and are therefore added for you. Please remove them from the -cov list" );
                            }

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

        logger.info( "The covariates being used here: " );
        logger.info( requestedCovariates );
    }
    
    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));
        // Only use data from non-dbsnp sites
        // Assume every mismatch at a non-dbsnp site is indicitive of poor quality
        if( dbsnp == null ) {
            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();
            SAMRecord read; // preallocate for use in for loop below
            int offset; // preallocate for use in for loop below
            for( int iii = 0; iii < reads.size(); iii++ ) {
                read = reads.get(iii);
                offset = offsets.get(iii);

                // Only use data from reads with mapping quality above given quality value and base quality greater than zero
                byte[] quals = read.getBaseQualities();
                if( read.getMappingQuality() >= MIN_MAPPING_QUALITY && quals[offset] > 0)
                {                    
                    if( offset > 0 && offset < (read.getReadLength() - 1) ) { // skip first and last bases because they don't have a dinuc count
                        updateDataFromRead(read, offset, ref);
                    }
                }

            }

        }
        
        return 1;
    }

    private void updateDataFromRead(SAMRecord read, int offset, ReferenceContext ref) {

        ArrayList<Comparable<?>> key = new ArrayList<Comparable<?>>();
        Comparable<?> keyElement; // preallocate for use in for loop below
        boolean badKey = false;
        for( Covariate covariate : requestedCovariates ) {
            keyElement = covariate.getValue( read, offset, ref.getBases() );
            if( keyElement != null ) {
                key.add( keyElement );
            } else {
                badKey = true;
            }
        }

        RecalDatum datum = null;
        if( !badKey ) {
            datum = dataManager.data.get( key );
            if( datum == null ) { // key doesn't exist yet in the map so make a new bucket and add it
                datum = new RecalDatum();
                dataManager.data.put( key, datum );
            }
        }

        byte[] bases = read.getReadBases();

        char base = (char)bases[offset];
        char refBase = ref.getBase();
        if ( read.getReadNegativeStrandFlag() ) {
            refBase = BaseUtils.simpleComplement( refBase );
            base = BaseUtils.simpleComplement( base );
        }

        if( datum != null ) {
            datum.increment( base, refBase );
        }
    }
   
    public PrintStream reduceInit() {
        try {
            return new PrintStream( RECAL_FILE );
        } catch ( FileNotFoundException e ) {
            throw new RuntimeException( "Couldn't open output file: ", e );
        }
    }

    public PrintStream reduce( Integer value, PrintStream recalTableStream ) {
        return recalTableStream; // nothing to do here
    }

    public void onTraversalDone( PrintStream recalTableStream ) {
        out.print( "Writing raw recalibration data..." );
        for( Covariate cov : requestedCovariates ) {
            recalTableStream.println( "@!" + cov.getClass().getSimpleName() ); // The "@!" is a code for TableRecalibrationWalker to recognize this line as a Covariate class name
        }
        outputToCSV( recalTableStream );
        out.println( "...done!" );
    
        recalTableStream.close();
    }

    private void outputToCSV( PrintStream recalTableStream ) {
        for( Map.Entry<List<? extends Comparable<?>>, RecalDatum> entry : dataManager.data.entrySet() ) {
            //for( int iii = 0; iii < entry.getKey().size(); iii++ ) {
                //index = Integer.parseInt( (entry.getKey().get( iii )).toString());
                //recalTableStream.print( requestedCovariates.get( iii ).getValue( index ) + "," );
            //	recalTableStream.print( entry.getKey().get(iii) )
            //}
        	for( Comparable<?> comp : entry.getKey() ) {
        		recalTableStream.print( comp + "," );
        	}
            recalTableStream.println( entry.getValue().outputToCSV() );
        }
    }

}

