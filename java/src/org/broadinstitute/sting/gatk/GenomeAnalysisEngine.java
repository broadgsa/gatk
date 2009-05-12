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
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.io.File;

public class GenomeAnalysisEngine {

    // our instance of this genome analysis toolkit; it's used by other classes to extract the traversal engine
    // TODO: public static without final tends to indicate we're thinking about this the wrong way
    public static GenomeAnalysisEngine instance;

    // our traversal engine
    private TraversalEngine engine = null;

    // the level of debugging we're using
    public boolean DEBUGGING = false;

    // our argument collection
    private final GATKArgumentCollection argCollection;

    /** Collection of output streams used by the walker. */
    private OutputTracker outputTracker = null;

    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(GenomeAnalysisEngine.class);

    /**
     * our constructor, where all the work is done
     * <p/>
     * legacy traversal types are sent to legacyTraversal function; as we move more of the traversals to the
     * new MicroScheduler class we'll be able to delete that function.
     *
     * @param args      the argument collection, where we get all our setup information from
     * @param my_walker the walker we're running with
     */
    public GenomeAnalysisEngine(GATKArgumentCollection args, Walker my_walker) {

        // validate our parameters
        if (args == null || my_walker == null) {
            throw new StingException("Neither the GATKArgumentCollection or the Walker passed to GenomeAnalysisEngine can be null.");
        }

        // save our argument parameter
        this.argCollection = args;

        // make sure our instance variable points to this analysis engine
        instance = this;

        // our reference ordered data collection
        List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods = new ArrayList<ReferenceOrderedData<? extends ReferenceOrderedDatum>>();

        //
        // please don't use these in the future, use the new syntax <- if we're not using these please remove them
        //
        if (argCollection.DBSNPFile != null) bindConvenienceRods("dbSNP", "dbsnp", argCollection.DBSNPFile);
        if (argCollection.HAPMAPFile != null)
            bindConvenienceRods("hapmap", "HapMapAlleleFrequencies", argCollection.HAPMAPFile);
        if (argCollection.HAPMAPChipFile != null)
            bindConvenienceRods("hapmap-chip", "GFF", argCollection.HAPMAPChipFile);

        // parse out the rod bindings
        ReferenceOrderedData.parseBindings(logger, argCollection.RODBindings, rods);

        // create the output streams
        initializeOutputStreams( my_walker );

        // our microscheduler, which is in charge of running everything
        MicroScheduler microScheduler = null;

        // if we're a read or a locus walker, we use the new system.  Right now we have complicated
        // branching based on the input data, but this should disapear when all the traversals are switched over
        if ((my_walker instanceof LocusWalker && argCollection.walkAllLoci && !(argCollection.samFiles == null || argCollection.samFiles.size() == 0)) ||
                my_walker instanceof ReadWalker) {
            microScheduler = createMicroscheduler(my_walker, rods);
        } else { // we have an old style traversal, once we're done return
            legacyTraversal(my_walker, rods);
            return;
        }

        // Prepare the sort ordering w.r.t. the sequence dictionary
        if (argCollection.referenceFile != null) {
            final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(argCollection.referenceFile);
            GenomeLoc.setupRefContigOrdering(refFile);
        }

        // Determine the validation stringency.  Default to ValidationStringency.STRICT.
        ValidationStringency strictness = getValidationStringency();

        logger.info("Strictness is " + strictness);

        // perform validation steps that are common to all the engines
        genericEngineSetup(strictness);

        // parse out any genomic location they've provided
        List<GenomeLoc> locs = setupIntervalRegion();

        // excute the microscheduler
        microScheduler.execute(my_walker, locs);
    }


    /**
     * this is to accomdate the older style traversals, that haven't been converted over to the new system.  Putting them
     * into their own function allows us to deviate in the two behaviors so the new style traversals aren't limited to what
     * the old style does.  As traversals are converted, this function should disappear.
     *
     * @param my_walker
     * @param rods
     */
    private void legacyTraversal(Walker my_walker, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        if (my_walker instanceof LocusWindowWalker) {
            this.engine = new TraverseByLocusWindows(argCollection.samFiles, argCollection.referenceFile, rods);
        } else if (my_walker instanceof LocusWalker) {
            if (argCollection.referenceFile == null)
                Utils.scareUser(String.format("Locus-based traversals require a reference file but none was given"));
            if (argCollection.samFiles == null || argCollection.samFiles.size() == 0) {
                if (((LocusWalker) my_walker).requiresReads())
                    Utils.scareUser(String.format("Analysis %s requires reads, but none were given", argCollection.analysisName));
                this.engine = new TraverseByReference(null, argCollection.referenceFile, rods);
            } else {
                if (((LocusWalker) my_walker).cannotHandleReads())
                    Utils.scareUser(String.format("Analysis %s doesn't support SAM/BAM reads, but a read file %s was provided", argCollection.analysisName, argCollection.samFiles));
                this.engine = new TraverseByLoci(argCollection.samFiles, argCollection.referenceFile, rods);
            }
        } else if (my_walker instanceof DuplicateWalker) {
            // we're a duplicate walker
            this.engine = new TraverseDuplicates(argCollection.samFiles, argCollection.referenceFile, rods);
        } else {
            throw new RuntimeException("Unexpected walker type: " + my_walker);
        }

        // Prepare the sort ordering w.r.t. the sequence dictionary
        if (argCollection.referenceFile != null) {
            final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(argCollection.referenceFile);
            GenomeLoc.setupRefContigOrdering(refFile);
        }

        // Determine the validation stringency.  Default to ValidationStringency.STRICT.
        ValidationStringency strictness = getValidationStringency();

        logger.info("Strictness is " + strictness);
        genericEngineSetup(strictness);


        engine.traverse(my_walker);

    }

    /**
     * setup a microscheduler
     *
     * @param my_walker our walker of type LocusWalker
     * @param rods      the reference order data
     * @return a new microscheduler
     */
    private MicroScheduler createMicroscheduler(Walker my_walker, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        // the mircoscheduler to return
        MicroScheduler microScheduler = null;

        // we need to verify different parameter based on the walker type
        if (my_walker instanceof LocusWalker) {
            // some warnings
            if (argCollection.referenceFile == null)
                Utils.scareUser(String.format("Locus-based traversals require a reference file but none was given"));
            if (((LocusWalker) my_walker).cannotHandleReads())
                Utils.scareUser(String.format("Analysis %s doesn't support SAM/BAM reads, but a read file %s was provided", argCollection.analysisName, argCollection.samFiles));

            // create the MicroScheduler
            microScheduler = MicroScheduler.create(my_walker, argCollection.samFiles, argCollection.referenceFile, rods, argCollection.numberOfThreads);
            engine = microScheduler.getTraversalEngine();
        }
        else if (my_walker instanceof ReadWalker)
        {
            if (argCollection.referenceFile == null)
                Utils.scareUser(String.format("Locus-based traversals require a reference file but none was given"));
            microScheduler = MicroScheduler.create(my_walker, argCollection.samFiles, argCollection.referenceFile, rods, argCollection.numberOfThreads);
            engine = microScheduler.getTraversalEngine();
        }

        return microScheduler;
    }


    /**
     * commands that get executed for each engine, regardless of the type
     *
     * @param strictness our current strictness level
     */
    private void genericEngineSetup(ValidationStringency strictness) {
        engine.setStrictness(strictness);

        engine.setMaxReads(Integer.parseInt(argCollection.maximumReads));

        // we default interval files over the genome region strin
        if (argCollection.intervals != null) {
            engine.setLocation(setupIntervalRegion());
        }
        // hmm...
        if (argCollection.maximumReadSorts != null) {
            engine.setSortOnFly(Integer.parseInt(argCollection.maximumReadSorts));
        }

        if (argCollection.downsampleFraction != null) {
            engine.setDownsampleByFraction(Double.parseDouble(argCollection.downsampleFraction));
        }

        if (argCollection.downsampleCoverage != null) {
            engine.setDownsampleByCoverage(Integer.parseInt(argCollection.downsampleCoverage));
        }

        engine.setSafetyChecking(!argCollection.unsafe);
        engine.setThreadedIO(argCollection.enabledThreadedIO);
        engine.setWalkOverAllSites(argCollection.walkAllLoci);
        engine.initialize();
    }


    /**
     * setup the interval regions, from either the interval file of the genome region string
     *
     * @return a list of genomeLoc representing the interval file
     */
    private List<GenomeLoc> setupIntervalRegion() {
        List<GenomeLoc> locs;
        if( new File(argCollection.intervals).exists() ) {
            logger.info("Intervals argument specifies a file.  Loading intervals from file.");
            return GenomeLoc.IntervalFileToList(argCollection.intervals);
        }
        else {
            logger.info("Intervals argument does not specify a file.  Trying to parse it as a simple string.");
            return GenomeLoc.parseGenomeLocs(argCollection.intervals);
        }
    }

    /**
     * Default to ValidationStringency.STRICT.
     *
     * @return the validation stringency
     */
    private ValidationStringency getValidationStringency() {
        ValidationStringency strictness;
        try {
            strictness = Enum.valueOf(ValidationStringency.class, argCollection.strictnessLevel);
        }
        catch (IllegalArgumentException ex) {
            strictness = ValidationStringency.STRICT;
        }
        return strictness;
    }

    /**
     * Convenience function that binds RODs using the old-style command line parser to the new style list for
     * a uniform processing.
     *
     * @param name
     * @param type
     * @param file
     */
    private void bindConvenienceRods(final String name, final String type, final String file) {
        argCollection.RODBindings.add(Utils.join(",", new String[]{name, type, file}));
    }


    /** Initialize the output streams as specified by the user. */
    private void initializeOutputStreams( Walker walker ) {
        outputTracker = (argCollection.outErrFileName != null) ? new OutputTracker(argCollection.outErrFileName, argCollection.outErrFileName)
                : new OutputTracker(argCollection.outFileName, argCollection.errFileName);
        walker.initializeOutputStreams(outputTracker);
    }

    /**
     * Gets the output tracker.  Tracks data available to a given walker.
     *
     * @return The output tracker.
     */
    public OutputTracker getOutputTracker() {
        return outputTracker;
    }


    public SAMFileReader getSamReader() {
        return this.engine.getSamReader();
    }

    public TraversalEngine getEngine() {
        return this.engine;
    }
}
