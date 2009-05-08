package org.broadinstitute.sting.gatk;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.executive.MicroScheduler;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.traversals.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class GenomeAnalysisTK extends CommandLineProgram {
    public static GenomeAnalysisTK Instance = null;

    // parameters and their defaults
    @Argument(fullName="input_file",shortName="I",doc="SAM or BAM file(s)",required=false)    
    protected List<File> INPUT_FILES = null;
    
    @Argument(fullName="maximum_reads",shortName="M",doc="Maximum number of reads to process before exiting",required=false)    
    protected String MAX_READS_ARG = "-1";

    @Argument(fullName="validation_strictness",shortName="S",doc="How strict should we be with validation (LENIENT|SILENT|STRICT)",required=false)    
    protected String STRICTNESS_ARG = "strict";
    
    @Argument(fullName="reference_sequence", shortName="R",doc="Reference sequence file",required=false)    
    protected File REF_FILE_ARG = null;

    @Argument(fullName="genome_region", shortName="L", doc="Genome region to operation on: from chr:start-end",required=false)    
    protected String REGION_STR = null;

    @Argument(fullName="analysis_type", shortName="T", doc="Type of analysis to run")    
    protected String Analysis_Name = null;

    @Argument(fullName="DBSNP",shortName="D",doc="DBSNP file",required=false)    
    protected String DBSNP_FILE = null;

    @Argument(fullName="hapmap",shortName="H",doc="Hapmap file",required=false)    
    protected String HAPMAP_FILE = null;

    @Argument(fullName="hapmap_chip",shortName="hc",doc="Hapmap chip file",required=false)
    protected String HAPMAP_CHIP_FILE = null;

    @Argument(fullName="threaded_IO",shortName="P",doc="If set, enables threaded I/O operations",required=false)    
    protected Boolean ENABLED_THREADED_IO = false;

    @Argument(fullName="unsafe",shortName="U",doc="If set, enables unsafe operations, nothing will be checked at runtime.",required=false)
    protected Boolean UNSAFE = false;

    @Argument(fullName="sort_on_the_fly",shortName="sort",doc="Maximum number of reads to sort on the fly",required=false)    
    protected String MAX_ON_FLY_SORTS = null;

    @Argument(fullName="downsample_to_fraction",shortName="dfrac",doc="Fraction [0.0-1.0] of reads to downsample to",required=false)
    protected String DOWNSAMPLE_FRACTION = null;

    @Argument(fullName="downsample_to_coverage",shortName="dcov",doc="Coverage [integer] to downsample to",required=false)
    protected String DOWNSAMPLE_COVERAGE = null;

    @Argument(fullName="intervals_file",shortName="V",doc="File containing list of genomic intervals to operate on. line := <contig> <start> <end>",required=false)    
    protected String INTERVALS_FILE = null;

    // our walker manager
    private WalkerManager walkerManager = null;

    public String pluginPathName = null;
    private TraversalEngine engine = null;
    public boolean DEBUGGING = false;

    @Argument(fullName="all_loci",shortName="A",doc="Should we process all loci, not just those covered by reads",required=false)    
    protected Boolean WALK_ALL_LOCI = false;

    @Argument(fullName="disablethreading",shortName="dt",doc="Disable experimental threading support.",required=false)    
    protected Boolean DISABLE_THREADING = false;

    /**
     * An output file presented to the walker.
     */
    @Argument(fullName="out",shortName="o",doc="An output file presented to the walker.  Will overwrite contents if file exists.",required=false)
    protected String outFileName = null;

    /**
     * An error output file presented to the walker.
     */
    @Argument(fullName="err",shortName="e",doc="An error output file presented to the walker.  Will overwrite contents if file exists.",required=false)
    protected String errFileName = null;

    /**
     * A joint file for both 'normal' and error output presented to the walker.
     */
    @Argument(fullName="outerr",shortName="oe",doc="A joint file for 'normal' and error output presented to the walker.  Will overwrite contents if file exists.",required=false)
    protected String outErrFileName = null;

    /**
     * How many threads should be allocated to this analysis.
     */
    @Argument(fullName="numthreads",shortName="nt",doc="How many threads should be allocated to running this analysis.",required=false)    
    protected int numThreads = 1;

    @Argument(fullName="rodBind",shortName="B",doc="",required=false)
    protected ArrayList<String> ROD_BINDINGS = new ArrayList<String>();    

    /**
     * Collection of output streams used by the walker.
     */
    private OutputTracker outputTracker = null;

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(GenomeAnalysisTK.class);



    /**
     * GATK can add arguments dynamically based on analysis type.
     * @return true
     */
    @Override
    protected boolean canAddArgumentsDynamically() { return true; }

    /**
     * GATK provides the walker as an argument source.  As a side-effect, initializes the walker variable.
     * @return List of walkers to load dynamically.
     */
    @Override
    protected Class[] getArgumentSources() {
        if( Analysis_Name == null )
            throw new IllegalArgumentException("Must provide analysis name");

        walkerManager = new WalkerManager( pluginPathName );

        if( !walkerManager.doesWalkerExist(Analysis_Name) )
            throw new IllegalArgumentException("Invalid analysis name");

        return new Class[] { walkerManager.getWalkerClassByName(Analysis_Name) };
    }

    @Override
    protected String getArgumentSourceName( Class argumentSource ) {
        return WalkerManager.getWalkerName( (Class<Walker>)argumentSource );
    }

    /**
     * Required main method implementation.
     */
    public static void main(String[] argv) {
        try {
            Instance = new GenomeAnalysisTK();
            start(Instance, argv);
        } catch ( Exception e ) {
            exitSystemWithError(e);
        }
    }

    /**
     * Convenience function that binds RODs using the old-style command line parser to the new style list for
     * a uniform processing.
     *
     * @param name
     * @param type
     * @param file
     */
    private void bindConvenienceRods(final String name, final String type, final String file )
    {
        ROD_BINDINGS.add( Utils.join(",", new String[] { name, type, file }) );
    }

    private static void printExitSystemMsg(final String msg) {
        System.out.printf("------------------------------------------------------------------------------------------%n");
        System.out.printf("An error has occurred%n");
        System.out.printf("It's possible (maybe even likely) that this is an input error on your part%n");
        System.out.printf("But if it's a bug or something that should work, please report this to gsadevelopers@broad.mit.edu%n");
        System.out.printf("%n");
        System.out.printf("%s%n", msg);
    }

    public static void exitSystemWithError(final String msg) {
        printExitSystemMsg(msg);
        System.exit(1);
    }

    public static void exitSystemWithError(final String msg, Exception e ) {
        e.printStackTrace();
        printExitSystemMsg(msg);
        System.exit(1);
    }

    public static void exitSystemWithError(Exception e ) {
        exitSystemWithError(e.getMessage(), e);
    }

    protected int execute() {
        final boolean TEST_ROD = false;
        List<ReferenceOrderedData<? extends ReferenceOrderedDatum> > rods = new ArrayList<ReferenceOrderedData<? extends ReferenceOrderedDatum> >();

        //
        // please don't use these in the future, use the new syntax
        //
        if ( DBSNP_FILE != null )               bindConvenienceRods("dbSNP", "dbsnp", DBSNP_FILE);
        if ( HAPMAP_FILE != null )              bindConvenienceRods("hapmap", "HapMapAlleleFrequencies", HAPMAP_FILE);
        if ( HAPMAP_CHIP_FILE != null )         bindConvenienceRods("hapmap-chip", "GFF", HAPMAP_CHIP_FILE);

        ReferenceOrderedData.parseBindings(logger, ROD_BINDINGS, rods);

        initializeOutputStreams();

        Walker<?,?> my_walker = null;
        try {
            my_walker = walkerManager.createWalkerByName( Analysis_Name );
            loadArgumentsIntoObject( my_walker );
        }
        catch( InstantiationException ex ) {
            throw new RuntimeException( "Unable to instantiate walker.", ex );
        }
        catch( IllegalAccessException ex ) {
            throw new RuntimeException( "Unable to access walker", ex );
        }

        MicroScheduler microScheduler = null;

        // Get the walker specified
        if ( my_walker instanceof LocusWalker ) {
            LocusWalker<?, ?> walker = (LocusWalker<?, ?>) my_walker;

            if ( REF_FILE_ARG == null )
                Utils.scareUser(String.format("Locus-based traversals require a reference file but none was given"));

            if ( INPUT_FILES == null || INPUT_FILES.size() == 0 ) {
                if ( walker.requiresReads() )
                    Utils.scareUser(String.format("Analysis %s requires reads, but none were given", Analysis_Name));
                this.engine = new TraverseByReference(null, REF_FILE_ARG, rods);
            } else {
                if ( walker.cannotHandleReads() )
                    Utils.scareUser(String.format("Analysis %s doesn't support SAM/BAM reads, but a read file %s was provided", Analysis_Name, INPUT_FILES));

                if ( WALK_ALL_LOCI ) {
                    microScheduler = MicroScheduler.create( walker, INPUT_FILES, REF_FILE_ARG, rods, numThreads );
                    this.engine = microScheduler.getTraversalEngine();
                }
                else
                	this.engine = new TraverseByLoci(INPUT_FILES, REF_FILE_ARG, rods);
            }
        } else if ( my_walker instanceof LocusWindowWalker ) {
            this.engine = new TraverseByLocusWindows(INPUT_FILES, REF_FILE_ARG, rods);
        } else if ( my_walker instanceof ReadWalker) {
            // we're a read walker
            this.engine = new TraverseByReads(INPUT_FILES, REF_FILE_ARG, rods);
        } else if ( my_walker instanceof DuplicateWalker) {
             // we're a duplicate walker
             this.engine = new TraverseDuplicates(INPUT_FILES, REF_FILE_ARG, rods);
        } else {
            throw new RuntimeException("Unexpected walker type: " + my_walker);
        }

        // Prepare the sort ordering w.r.t. the sequence dictionary
        if (REF_FILE_ARG != null) {
            final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REF_FILE_ARG);
            GenomeLoc.setupRefContigOrdering(refFile);
        }

        // Determine the validation stringency.  Default to ValidationStringency.STRICT.
        ValidationStringency strictness;
        try {
            strictness = Enum.valueOf(ValidationStringency.class, STRICTNESS_ARG);
        }
        catch( IllegalArgumentException ex ) {
            strictness = ValidationStringency.STRICT;
        }

        logger.info("Strictness is " + strictness);
        engine.setStrictness(strictness);

        engine.setMaxReads(Integer.parseInt(MAX_READS_ARG));

        if (REGION_STR != null) {
            engine.setLocation(REGION_STR);
        }

        if (INTERVALS_FILE != null) {
            engine.setLocationFromFile(INTERVALS_FILE);
        }

        if (MAX_ON_FLY_SORTS != null) {
            engine.setSortOnFly(Integer.parseInt(MAX_ON_FLY_SORTS));
        }

        if (DOWNSAMPLE_FRACTION != null) {
            engine.setDownsampleByFraction(Double.parseDouble(DOWNSAMPLE_FRACTION));
        }

        if (DOWNSAMPLE_COVERAGE != null) {
            engine.setDownsampleByCoverage(Integer.parseInt(DOWNSAMPLE_COVERAGE));
        }

        engine.setSafetyChecking(!UNSAFE);
        engine.setThreadedIO(ENABLED_THREADED_IO);
        engine.setWalkOverAllSites(WALK_ALL_LOCI);
        engine.initialize();

        if( microScheduler != null ) {
            List<GenomeLoc> locs = null;
            if (INTERVALS_FILE != null)
                locs = GenomeLoc.IntervalFileToList(INTERVALS_FILE);
            else
                locs = GenomeLoc.parseGenomeLocs( REGION_STR );
            microScheduler.execute( my_walker, locs );
        }
        else
            engine.traverse(my_walker);

        return 0;
    }

    /**
     * Initialize the output streams as specified by the user.
     */
    private void initializeOutputStreams() {
        outputTracker = (outErrFileName != null) ? new OutputTracker( outErrFileName, outErrFileName )
                                                 : new OutputTracker( outFileName, errFileName );
    }

    /**
     * Gets the output tracker.  Tracks data available to a given walker.
     * @return The output tracker.
     */
    public OutputTracker getOutputTracker() {
        return outputTracker;
    }


    public SAMFileReader getSamReader() { return this.engine.getSamReader(); }
    public TraversalEngine getEngine() { return this.engine; }
}
