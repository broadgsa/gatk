/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine;

import com.google.java.contract.Ensures;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.engine.filters.DisableableReadFilter;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.engine.datasources.reads.*;
import org.broadinstitute.gatk.engine.datasources.reference.ReferenceDataSource;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.gatk.engine.executive.MicroScheduler;
import org.broadinstitute.gatk.engine.filters.FilterManager;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.filters.ReadGroupBlackListFilter;
import org.broadinstitute.gatk.engine.io.OutputTracker;
import org.broadinstitute.gatk.engine.io.stubs.Stub;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.iterators.ReadTransformersMode;
import org.broadinstitute.gatk.engine.phonehome.GATKRunReport;
import org.broadinstitute.gatk.utils.io.ReferenceBacked;
import org.broadinstitute.gatk.utils.refdata.tracks.IndexDictionaryUtils;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.engine.samples.SampleDB;
import org.broadinstitute.gatk.engine.samples.SampleDBBuilder;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.genotyper.IndexedSampleList;
import org.broadinstitute.gatk.utils.genotyper.SampleList;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import org.broadinstitute.gatk.utils.progressmeter.ProgressMeter;
import org.broadinstitute.gatk.engine.recalibration.BQSRArgumentSet;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.utils.text.XReadLines;
import org.broadinstitute.gatk.utils.threading.ThreadEfficiencyMonitor;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.concurrent.TimeUnit;

import static org.broadinstitute.gatk.utils.DeprecatedToolChecks.getWalkerDeprecationInfo;
import static org.broadinstitute.gatk.utils.DeprecatedToolChecks.isDeprecatedWalker;

/**
 * A GenomeAnalysisEngine that runs a specified walker.
 */
public class GenomeAnalysisEngine {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(GenomeAnalysisEngine.class);
    public static final long NO_RUNTIME_LIMIT = -1;

    /**
     * The GATK command-line argument parsing code.
     */
    private ParsingEngine parsingEngine;

    /**
     * The genomeLocParser can create and parse GenomeLocs.
     */
    private GenomeLocParser genomeLocParser;

    /**
     * Accessor for sharded read data.
     */
    private SAMDataSource readsDataSource = null;

    /**
     * Accessor for sharded reference data.
     */
    private ReferenceDataSource referenceDataSource = null;

    /**
     * Accessor for sample metadata
     */
    private SampleDB sampleDB = new SampleDB();

    /**
     * Accessor for sharded reference-ordered data.
     */
    private List<ReferenceOrderedDataSource> rodDataSources;

    // our argument collection
    private GATKArgumentCollection argCollection;

    /**
     * Collection of intervals used by the engine.
     */
    private GenomeLocSortedSet intervals = null;

    /**
     * Explicitly assign the interval set to use for this traversal (for unit testing purposes)
     * @param intervals set of intervals to use for this traversal
     */
    public void setIntervals( GenomeLocSortedSet intervals ) {
        this.intervals = intervals;
    }

    /**
     * Collection of inputs used by the engine.
     */
    private Map<ArgumentSource, Object> inputs = new HashMap<ArgumentSource, Object>();

    /**
     * Collection of outputs used by the engine.
     */
    private Collection<Stub<?>> outputs = new ArrayList<Stub<?>>();

    /**
     * Collection of the filters applied to the input data.
     */
    private Collection<ReadFilter> filters;

    /**
     * Collection of the read transformers applied to the reads
     */
    private List<ReadTransformer> readTransformers;

    /**
     * Controls the allocation of threads between CPU vs IO.
     */
    private ThreadAllocation threadAllocation;

    private ReadMetrics cumulativeMetrics = null;

    /**
     * A currently hacky unique name for this GATK instance
     */
    private String myName = "GATK_" + Math.abs(Utils.getRandomGenerator().nextInt());

    /**
     * our walker manager
     */
    private final WalkerManager walkerManager = new WalkerManager();

    private Walker<?, ?> walker;

    public void setWalker(Walker<?, ?> walker) {
        this.walker = walker;
    }

    /**
     * The short name of the current GATK walker as a string
     * @return a non-null String
     */
    public String getWalkerName() {
        return getWalkerName(walker.getClass());
    }

    /**
     * A processed collection of SAM reader identifiers.
     */
    private Collection<SAMReaderID> samReaderIDs = Collections.emptyList();

    /**
     * Set the SAM/BAM files over which to traverse.
     * @param samReaderIDs Collection of ids to use during this traversal.
     */
    public void setSAMFileIDs(Collection<SAMReaderID> samReaderIDs) {
        this.samReaderIDs = samReaderIDs;
    }

    /**
     * Collection of reference metadata files over which to traverse.
     */
    private Collection<RMDTriplet> referenceMetaDataFiles;

    /**
     * The threading efficiency monitor we use in the GATK to monitor our efficiency.
     *
     * May be null if one isn't active, or hasn't be initialized yet
     */
    private ThreadEfficiencyMonitor threadEfficiencyMonitor = null;

    /**
     * The global progress meter we are using to track our progress through the genome
     */
    private ProgressMeter progressMeter = null;

    /**
     * Set the reference metadata files to use for this traversal.
     * @param referenceMetaDataFiles Collection of files and descriptors over which to traverse.
     */
    public void setReferenceMetaDataFiles(Collection<RMDTriplet> referenceMetaDataFiles) {
        this.referenceMetaDataFiles = referenceMetaDataFiles;
    }

    /**
     * The maximum runtime of this engine, in nanoseconds, set during engine initialization
     * from the GATKArgumentCollection command line value
     */
    private long runtimeLimitInNanoseconds = -1;

    /**
     *  Base Quality Score Recalibration helper object
     */
    private BQSRArgumentSet bqsrArgumentSet = null;
    public BQSRArgumentSet getBQSRArgumentSet() { return bqsrArgumentSet; }
    public boolean hasBQSRArgumentSet() { return bqsrArgumentSet != null; }
    public void setBaseRecalibration(final GATKArgumentCollection args) {
        bqsrArgumentSet = new BQSRArgumentSet(args);
    }

    /**
     * Actually run the GATK with the specified walker.
     *
     * @return the value of this traversal.
     */
    public Object execute() {
        // first thing is to make sure the AWS keys can be decrypted
        GATKRunReport.checkAWSAreValid();

        //HeapSizeMonitor monitor = new HeapSizeMonitor();
        //monitor.start();
        setStartTime(new java.util.Date());

        final GATKArgumentCollection args = this.getArguments();

        // validate our parameters
        if (args == null) {
            throw new ReviewedGATKException("The GATKArgumentCollection passed to GenomeAnalysisEngine cannot be null.");
        }

        // validate our parameters              
        if (this.walker == null)
            throw new ReviewedGATKException("The walker passed to GenomeAnalysisEngine cannot be null.");

        // check that active region walkers do not use the downsampling to coverage argument
        checkDownSamplingToCoverage();
        
        if (args.nonDeterministicRandomSeed)
            Utils.resetRandomGenerator(System.currentTimeMillis());

        // if the use specified an input BQSR recalibration table then enable on the fly recalibration
        if (args.BQSR_RECAL_FILE != null)
            setBaseRecalibration(args);

        // setup the runtime limits
        setupRuntimeLimits(args);

        // Determine how the threads should be divided between CPU vs. IO.
        determineThreadAllocation();

        // Prepare the data for traversal.
        initializeDataSources();

        // initialize and validate the interval list
        initializeIntervals();
        validateSuppliedIntervals();

        // check to make sure that all sequence dictionaries are compatible with the reference's sequence dictionary
        validateDataSourcesAgainstReference(readsDataSource, referenceDataSource.getReference(), rodDataSources);

        // initialize sampleDB
        initializeSampleDB();

        // our microscheduler, which is in charge of running everything
        MicroScheduler microScheduler = createMicroscheduler();
        threadEfficiencyMonitor = microScheduler.getThreadEfficiencyMonitor();

        // create temp directories as necessary
        initializeTempDirectory();

        // create the output streams
        initializeOutputStreams(microScheduler.getOutputTracker());

        // Initializing the shard iterator / BAM schedule might take some time, so let the user know vaguely what's going on
        logger.info("Preparing for traversal" +
                    (readsDataSource.getReaderIDs().size() > 0 ? String.format(" over %d BAM files", readsDataSource.getReaderIDs().size()) : ""));
        Iterable<Shard> shardStrategy = getShardStrategy(readsDataSource,microScheduler.getReference(),intervals);
        logger.info("Done preparing for traversal");

        // execute the microscheduler, storing the results
        return microScheduler.execute(this.walker, shardStrategy);

        //monitor.stop();
        //logger.info(String.format("Maximum heap size consumed: %d",monitor.getMaxMemoryUsed()));

        //return result;
    }

    /**
     * Retrieves an instance of the walker based on the walker name.
     *
     * @param walkerName Name of the walker.  Must not be null.  If the walker cannot be instantiated, an exception will be thrown.
     * @return An instance of the walker.
     */
    public Walker<?, ?> getWalkerByName(String walkerName) {
        try {
            return walkerManager.createByName(walkerName);
        } catch ( UserException e ) {
            if ( isDeprecatedWalker(walkerName) ) {
                e = new UserException.DeprecatedWalker(walkerName, getWalkerDeprecationInfo(walkerName));
            }
            throw e;
        }
    }

    /**
     * Gets the name of a given walker type.
     * @param walkerType Type of walker.
     * @return Name of the walker.
     */
    public String getWalkerName(Class<? extends Walker> walkerType) {
        return walkerManager.getName(walkerType);
    }

    public String getName() {
        return myName;
    }

    /**
     * Gets a list of the filters to associate with the given walker.  Will NOT initialize the engine with this filters;
     * the caller must handle that directly.
     * @return A collection of available filters.
     */
    public Collection<ReadFilter> createFilters() {
        final List<ReadFilter> filters = new LinkedList<>();

        // First add the user requested filters
        if (this.getArguments().readGroupBlackList != null && !this.getArguments().readGroupBlackList.isEmpty())
            filters.add(new ReadGroupBlackListFilter(this.getArguments().readGroupBlackList));
        for(final String filterName: this.getArguments().readFilters)
            filters.add(this.getFilterManager().createByName(filterName));

        // now add the walker default filters.  This ordering is critical important if
        // users need to apply filters that fix up reads that would be removed by default walker filters
        filters.addAll(WalkerManager.getReadFilters(walker,this.getFilterManager()));

        // disable user-specified read filters, if allowed
        for(final String filterName: this.getArguments().disabledReadFilters) {
            ReadFilter filterToDisable = this.getFilterManager().createByName(filterName);
            if (! (filterToDisable instanceof DisableableReadFilter))
                throw new IllegalStateException(filterToDisable + " cannot be disabled");

            // so we're not trying to modify the list we're iterating over
            List<ReadFilter> filtersCopy = new ArrayList<>(filters);
            for (ReadFilter filter : filtersCopy) {
                if (filter.getClass() == filterToDisable.getClass())
                    filters.remove(filter);
            }
        }

        return Collections.unmodifiableList(filters);
    }

    /**
     * Returns a list of active, initialized read transformers
     *
     * @param walker the walker we need to apply read transformers too
     */
    public void initializeReadTransformers(final Walker walker) {
        // keep a list of the active read transformers sorted based on priority ordering
        List<ReadTransformer> activeTransformers = new ArrayList<ReadTransformer>();

        final ReadTransformersMode overrideMode = WalkerManager.getWalkerAnnotation(walker, ReadTransformersMode.class);
        final ReadTransformer.ApplicationTime overrideTime = overrideMode != null ? overrideMode.ApplicationTime() : null;

        final PluginManager<ReadTransformer> pluginManager = new PluginManager<ReadTransformer>(ReadTransformer.class);

        for ( final ReadTransformer transformer : pluginManager.createAllTypes() ) {
            transformer.initialize(overrideTime, this, walker);
            if ( transformer.enabled() )
                activeTransformers.add(transformer);
        }

        setReadTransformers(activeTransformers);
    }

    public List<ReadTransformer> getReadTransformers() {
        return readTransformers;
    }

    /*
     * Sanity checks that incompatible read transformers are not active together (and throws an exception if they are).
     *
     * @param readTransformers   the active read transformers
     */
    protected void checkActiveReadTransformers(final List<ReadTransformer> readTransformers) {
        if ( readTransformers == null )
            throw new IllegalArgumentException("read transformers cannot be null");

        ReadTransformer sawMustBeFirst = null;
        ReadTransformer sawMustBeLast  = null;

        for ( final ReadTransformer r : readTransformers ) {
            if ( r.getOrderingConstraint() == ReadTransformer.OrderingConstraint.MUST_BE_FIRST ) {
                if ( sawMustBeFirst != null )
                    throw new UserException.IncompatibleReadFiltersException(sawMustBeFirst.toString(), r.toString());
                sawMustBeFirst = r;
            } else if ( r.getOrderingConstraint() == ReadTransformer.OrderingConstraint.MUST_BE_LAST ) {
                if ( sawMustBeLast != null )
                    throw new UserException.IncompatibleReadFiltersException(sawMustBeLast.toString(), r.toString());
                sawMustBeLast = r;
            }
        }
    }

    protected void setReadTransformers(final List<ReadTransformer> readTransformers) {
        if ( readTransformers == null )
            throw new ReviewedGATKException("read transformers cannot be null");

        // sort them in priority order
        Collections.sort(readTransformers, new ReadTransformer.ReadTransformerComparator());

        // make sure we don't have an invalid set of active read transformers
        checkActiveReadTransformers(readTransformers);

        this.readTransformers = readTransformers;
    }

    /**
     * Parse out the thread allocation from the given command-line argument.
     */
    private void determineThreadAllocation() {
        if ( argCollection.numberOfDataThreads < 1 ) throw new UserException.BadArgumentValue("num_threads", "cannot be less than 1, but saw " + argCollection.numberOfDataThreads);
        if ( argCollection.numberOfCPUThreadsPerDataThread < 1 ) throw new UserException.BadArgumentValue("num_cpu_threads", "cannot be less than 1, but saw " + argCollection.numberOfCPUThreadsPerDataThread);
        if ( argCollection.numberOfIOThreads < 0 ) throw new UserException.BadArgumentValue("num_io_threads", "cannot be less than 0, but saw " + argCollection.numberOfIOThreads);

        this.threadAllocation = new ThreadAllocation(argCollection.numberOfDataThreads,
                argCollection.numberOfCPUThreadsPerDataThread,
                argCollection.numberOfIOThreads,
                argCollection.monitorThreadEfficiency);
    }

    public int getTotalNumberOfThreads() {
        return this.threadAllocation == null ? 1 : threadAllocation.getTotalNumThreads();
    }

    /**
     * Allow subclasses and others within this package direct access to the walker manager.
     * @return The walker manager used by this package.
     */
    protected WalkerManager getWalkerManager() {
        return walkerManager;
    }
    
    /**
     * setup a microscheduler
     *
     * @return a new microscheduler
     */
    private MicroScheduler createMicroscheduler() {
        // Temporarily require all walkers to have a reference, even if that reference is not conceptually necessary.
        if ((walker instanceof ReadWalker || walker instanceof DuplicateWalker || walker instanceof ReadPairWalker) &&
                this.getArguments().referenceFile == null) {
            throw new UserException.CommandLineException("Read-based traversals require a reference file but none was given");
        }

        return MicroScheduler.create(this,walker,this.getReadsDataSource(),this.getReferenceDataSource().getReference(),this.getRodDataSources(),threadAllocation);
    }

    protected DownsamplingMethod getDownsamplingMethod() {
        GATKArgumentCollection argCollection = this.getArguments();

        DownsamplingMethod commandLineMethod = argCollection.getDownsamplingMethod();
        DownsamplingMethod walkerMethod = WalkerManager.getDownsamplingMethod(walker);

        DownsamplingMethod method = commandLineMethod != null ? commandLineMethod : walkerMethod;
        checkCompatibilityWithWalker(method, walker);
        return method;
    }

    private static void checkCompatibilityWithWalker( DownsamplingMethod method, Walker walker ) {
        // Refactored from DownsamplingMethod
        final DownsampleType type = method.type;
        final Integer toCoverage = method.toCoverage;
        final boolean isLocusTraversal = walker instanceof LocusWalker || walker instanceof ActiveRegionWalker;

        if ( isLocusTraversal && type == DownsampleType.ALL_READS && toCoverage != null ) {
            throw new UserException("Downsampling to coverage with the ALL_READS method for locus-based traversals (eg., LocusWalkers) is not currently supported (though it is supported for ReadWalkers).");
        }

        // For locus traversals, ensure that the dcov value (if present) is not problematically low
        if ( isLocusTraversal && type != DownsampleType.NONE && toCoverage != null &&
                toCoverage < DownsamplingMethod.MINIMUM_SAFE_COVERAGE_TARGET_FOR_LOCUS_BASED_TRAVERSALS ) {
            throw new UserException(String.format("Locus-based traversals (ie., Locus and ActiveRegion walkers) require " +
                            "a minimum -dcov value of %d when downsampling to coverage. Values less " +
                            "than this can produce problematic downsampling artifacts while providing " +
                            "only insignificant improvements in memory usage in most cases.",
                    DownsamplingMethod.MINIMUM_SAFE_COVERAGE_TARGET_FOR_LOCUS_BASED_TRAVERSALS));
        }
    }

    protected void setDownsamplingMethod(DownsamplingMethod method) {
        argCollection.setDownsamplingMethod(method);
    }

    protected boolean includeReadsWithDeletionAtLoci() {
        return walker.includeReadsWithDeletionAtLoci();
    }

    /**
     * Verifies that the supplied set of reads files mesh with what the walker says it requires;
     * also makes sure that list of SAM files specified on the command line is not empty and contains
     * no duplicates.
     */
    protected void validateSuppliedReads() {
        GATKArgumentCollection arguments = this.getArguments();
        final Boolean samFilesArePresent = (arguments.samFiles != null && !arguments.samFiles.isEmpty());

        // Check what the walker says is required against what was provided on the command line.
        if (WalkerManager.isRequired(walker, DataSource.READS) && !samFilesArePresent)
            throw new ArgumentException("Walker requires reads but none were provided.");

        // Check what the walker says is allowed against what was provided on the command line.
        if (samFilesArePresent && !WalkerManager.isAllowed(walker, DataSource.READS))
            throw new ArgumentException("Walker does not allow reads but reads were provided.");

        //Make sure SAM list specified by the user (if necessary) is not empty
        if(WalkerManager.isRequired(walker, DataSource.READS) && samFilesArePresent && samReaderIDs.isEmpty() ) {
            throw new UserException("The list of input files does not contain any BAM files.");
        }

        // Make sure no SAM files were specified multiple times by the user.
        checkForDuplicateSamFiles();
    }

    /**
     * Check that active region walkers do not use the downsampling to coverage argument
     *
     * @throws UserException if an active region walker is using the -dcov or --downsample_to_coverage downsampling arguments
     */
    private void checkDownSamplingToCoverage() {
        if (argCollection.downsampleCoverage != null && walker instanceof ActiveRegionWalker) {
            throw new UserException.CommandLineException("Cannot use -dcov or --downsample_to_coverage for ActiveRegionWalkers, use another downsampling argument");
        }
    }

    /**
     * Checks whether there are SAM files that appear multiple times in the fully unpacked list of
     * SAM files (samReaderIDs). If there are, throws an ArgumentException listing the files in question.
     */
    protected void checkForDuplicateSamFiles() {
        Set<SAMReaderID> encounteredSamFiles = new HashSet<SAMReaderID>();
        Set<String> duplicateSamFiles = new LinkedHashSet<String>();

        for ( SAMReaderID samFile : samReaderIDs ) {
            if ( encounteredSamFiles.contains(samFile) ) {
                duplicateSamFiles.add(samFile.getSamFilePath());
            }
            else {
                encounteredSamFiles.add(samFile);
            }
        }

        if ( duplicateSamFiles.size() > 0 ) {
            throw new UserException("The following BAM files appear multiple times in the list of input files: " +
                                    duplicateSamFiles + " BAM files may be specified at most once.");
        }

    }

    /**
     * Verifies that the supplied reference file mesh with what the walker says it requires.
     */
    protected void validateSuppliedReference() {
        GATKArgumentCollection arguments = this.getArguments();
        // Check what the walker says is required against what was provided on the command line.
        // TODO: Temporarily disabling WalkerManager.isRequired check on the reference because the reference is always required.
        if (/*WalkerManager.isRequired(walker, DataSource.REFERENCE) &&*/ arguments.referenceFile == null)
            throw new ArgumentException("Walker requires a reference but none was provided.");

        // Check what the walker says is allowed against what was provided on the command line.
        if (arguments.referenceFile != null && !WalkerManager.isAllowed(walker, DataSource.REFERENCE))
            throw new ArgumentException("Walker does not allow a reference but one was provided.");
    }

    protected void validateSuppliedIntervals() {
        // Only read walkers support '-L unmapped' intervals.  Trap and validate any other instances of -L unmapped.
        if(!(walker instanceof ReadWalker)) {
            GenomeLocSortedSet intervals = getIntervals();
            if(intervals != null && getIntervals().contains(GenomeLoc.UNMAPPED))
                throw new ArgumentException("Interval list specifies unmapped region.  Only read walkers may include the unmapped region.");
        }

        // If intervals is non-null and empty at this point, it means that the list of intervals to process
        // was filtered down to an empty set (eg., the user specified something like -L chr1 -XL chr1). Since
        // this was very likely unintentional, the user should be informed of this. Note that this is different
        // from the case where intervals == null, which indicates that there were no interval arguments.
        if ( intervals != null && intervals.isEmpty() ) {
            logger.warn("The given combination of -L and -XL options results in an empty set.  No intervals to process.");
        }

        // TODO: add a check for ActiveRegion walkers to prevent users from passing an entire contig/chromosome
    }

    /**
     * Get the sharding strategy given a driving data source.
     *
     * @param readsDataSource readsDataSource
     * @param drivingDataSource Data on which to shard.
     * @param intervals intervals
     * @return the sharding strategy
     */
    protected Iterable<Shard> getShardStrategy(SAMDataSource readsDataSource, ReferenceSequenceFile drivingDataSource, GenomeLocSortedSet intervals) {
        ValidationExclusion exclusions = (readsDataSource != null ? readsDataSource.getReadsInfo().getValidationExclusionList() : null);
        DownsamplingMethod downsamplingMethod = readsDataSource != null ? readsDataSource.getReadsInfo().getDownsamplingMethod() : null;
        ReferenceDataSource referenceDataSource = this.getReferenceDataSource();

        // If reads are present, assume that accessing the reads is always the dominant factor and shard based on that supposition.
        if(!readsDataSource.isEmpty()) {
            if(!readsDataSource.hasIndex() && !exclusions.contains(ValidationExclusion.TYPE.ALLOW_UNINDEXED_BAM))
                throw new UserException.CommandLineException("Cannot process the provided BAM/CRAM file(s) because they were not indexed.  The GATK does offer limited processing of unindexed BAM/CRAMs in --unsafe mode, but this feature is unsupported -- use it at your own risk!");
            if(!readsDataSource.hasIndex() && intervals != null && !argCollection.allowIntervalsWithUnindexedBAM)
                throw new UserException.CommandLineException("Cannot perform interval processing when reads are present but no index is available.");

            if(walker instanceof LocusWalker) {
                if (readsDataSource.getSortOrder() != SAMFileHeader.SortOrder.coordinate)
                    throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.coordinate, "Locus walkers can only traverse coordinate-sorted data.  Please resort your input BAM file(s) or set the Sort Order tag in the header appropriately.");
                if(intervals == null)
                    return readsDataSource.createShardIteratorOverMappedReads(new LocusShardBalancer());
                else
                    return readsDataSource.createShardIteratorOverIntervals(intervals,new LocusShardBalancer());
            } 
            else if(walker instanceof ActiveRegionWalker) {
                if (readsDataSource.getSortOrder() != SAMFileHeader.SortOrder.coordinate)
                    throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.coordinate, "Active region walkers can only traverse coordinate-sorted data.  Please resort your input BAM file(s) or set the Sort Order tag in the header appropriately.");
                if(intervals == null)
                    return readsDataSource.createShardIteratorOverMappedReads(new ActiveRegionShardBalancer());
                else
                    return readsDataSource.createShardIteratorOverIntervals(((ActiveRegionWalker)walker).extendIntervals(intervals, this.genomeLocParser, this.getReferenceDataSource().getReference()), new ActiveRegionShardBalancer());
            } 
            else if(walker instanceof ReadWalker || walker instanceof ReadPairWalker || walker instanceof DuplicateWalker) {
                // Apply special validation to read pair walkers.
                if(walker instanceof ReadPairWalker) {
                    if(readsDataSource.getSortOrder() != SAMFileHeader.SortOrder.queryname)
                        throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.queryname, "Read pair walkers are exceptions in that they cannot be run on coordinate-sorted BAMs but instead require query name-sorted files.  You will need to resort your input BAM file in query name order to use this walker.");
                    if(intervals != null && !intervals.isEmpty())
                        throw new UserException.CommandLineException("Pairs traversal cannot be used in conjunction with intervals.");
                }

                if(intervals == null)
                    return readsDataSource.createShardIteratorOverAllReads(new ReadShardBalancer());
                else
                    return readsDataSource.createShardIteratorOverIntervals(intervals, new ReadShardBalancer());
            }
            else
                throw new ReviewedGATKException("Unable to determine walker type for walker " + walker.getClass().getName());
        }
        else {
            // TODO -- Determine what the ideal shard size should be here.  Matt suggested that a multiple of 16K might work well
            // TODO --  (because of how VCF indexes work), but my empirical experience has been simply that the larger the shard
            // TODO --  size the more efficient the traversal (at least for RODWalkers).  Keeping the previous values for now.  [EB]
            final int SHARD_SIZE = walker instanceof RodWalker ? 1000000 : 100000;
            if(intervals == null)
                return referenceDataSource.createShardsOverEntireReference(readsDataSource,genomeLocParser,SHARD_SIZE);
            else
                return referenceDataSource.createShardsOverIntervals(readsDataSource,intervals,SHARD_SIZE);
        }
    }

    protected boolean flashbackData() {
        return walker instanceof ReadWalker;
    }

    /**
     * Create the temp directory if it doesn't exist.
     */
    private void initializeTempDirectory() {
        File tempDir = new File(System.getProperty("java.io.tmpdir"));
        if (!tempDir.exists() && !tempDir.mkdirs())
            throw new UserException.BadTmpDir("Unable to create directory");
    }

    /**
     * Initialize the output streams as specified by the user.
     *
     * @param outputTracker the tracker supplying the initialization data.
     */
    private void initializeOutputStreams(final OutputTracker outputTracker) {
        for (final Map.Entry<ArgumentSource, Object> input : getInputs().entrySet()) {
            setReferenceFile(input.getValue());
            outputTracker.addInput(input.getKey(), input.getValue());
        }
        for (final Stub<?> stub : getOutputs()) {
            setReferenceFile(stub);
            stub.processArguments(argCollection);
            outputTracker.addOutput(stub);
        }

        outputTracker.prepareWalker(walker, getArguments().strictnessLevel);
    }

    private void setReferenceFile(final Object object) {
        if (object instanceof ReferenceBacked) {
            ((ReferenceBacked)object).setReferenceFile(argCollection.referenceFile);
        }
    }

    public ReferenceDataSource getReferenceDataSource() {
        return referenceDataSource;
    }

    public GenomeLocParser getGenomeLocParser() {
        return genomeLocParser;
    }

    /**
     * Manage lists of filters.
     */
    private final FilterManager filterManager = new FilterManager();

    private Date startTime = null; // the start time for execution

    public void setParser(ParsingEngine parsingEngine) {
        this.parsingEngine = parsingEngine;
    }

    /**
     * Explicitly set the GenomeLocParser, for unit testing.
     * @param genomeLocParser GenomeLocParser to use.
     */
    public void setGenomeLocParser(GenomeLocParser genomeLocParser) {
        this.genomeLocParser = genomeLocParser;
    }

    /**
     * Sets the start time when the execute() function was last called
     * @param startTime the start time when the execute() function was last called
     */
    protected void setStartTime(Date startTime) {
        this.startTime = startTime;
    }

    /**
     * @return the start time when the execute() function was last called
     */
    public Date getStartTime() {
        return startTime;
    }

    /**
     * Setup the intervals to be processed
     */
    protected void initializeIntervals() {
        intervals = IntervalUtils.parseIntervalArguments(this.referenceDataSource.getReference(), argCollection.intervalArguments);
    }

    /**
     * Add additional, externally managed IO streams for inputs.
     *
     * @param argumentSource Field into which to inject the value.
     * @param value          Instance to inject.
     */
    public void addInput(ArgumentSource argumentSource, Object value) {
        inputs.put(argumentSource, value);
    }

    /**
     * Add additional, externally managed IO streams for output.
     *
     * @param stub Instance to inject.
     */
    public void addOutput(Stub<?> stub) {
        outputs.add(stub);
    }

    /**
     * Returns the tag associated with a given command-line argument.
     * @param key Object for which to inspect the tag.
     * @return Tags object associated with the given key, or an empty Tag structure if none are present. 
     */
    public Tags getTags(Object key)  {
        return parsingEngine.getTags(key);
    }

    protected void initializeDataSources() {
        logger.info("Strictness is " + argCollection.strictnessLevel);

        validateSuppliedReference();
        setReferenceDataSource(argCollection.referenceFile);

        validateSuppliedReads();
        initializeReadTransformers(walker);

        final Map<String, String> sampleRenameMap = argCollection.sampleRenameMappingFile != null ?
                                                    loadSampleRenameMap(argCollection.sampleRenameMappingFile) :
                                                    null;

        readsDataSource = createReadsDataSource(argCollection,genomeLocParser,referenceDataSource.getReference(), sampleRenameMap);

        for (ReadFilter filter : filters)
            filter.initialize(this);

        // set the sequence dictionary of all of Tribble tracks to the sequence dictionary of our reference
        rodDataSources = getReferenceOrderedDataSources(referenceMetaDataFiles,referenceDataSource.getReference().getSequenceDictionary(),
                                                        genomeLocParser,argCollection.unsafe,sampleRenameMap);
    }

    /**
     * Purely for testing purposes.  Do not use unless you absolutely positively know what you are doing (or
     * need to absolutely positively kill everyone in the room)
     * @param dataSource
     */
    public void setReadsDataSource(final SAMDataSource dataSource) {
        this.readsDataSource = dataSource;
    }

    /**
     * Entry-point function to initialize the samples database from input data and pedigree arguments
     */
    private void initializeSampleDB() {
        SampleDBBuilder sampleDBBuilder = new SampleDBBuilder(this, argCollection.pedigreeValidationType);
        sampleDBBuilder.addSamplesFromSAMHeader(getSAMFileHeader());
        sampleDBBuilder.addSamplesFromSampleNames(SampleUtils.getUniqueSamplesFromRods(this));
        sampleDBBuilder.addSamplesFromPedigreeFiles(argCollection.pedigreeFiles);
        sampleDBBuilder.addSamplesFromPedigreeStrings(argCollection.pedigreeStrings);
        sampleDB = sampleDBBuilder.getFinalSampleDB();
    }

    /**
     * Gets a unique identifier for the reader sourcing this read.
     * @param read Read to examine.
     * @return A unique identifier for the source file of this read.  Exception if not found.
     */
    public SAMReaderID getReaderIDForRead(final SAMRecord read) {
        return getReadsDataSource().getReaderID(read);
    }

    /**
     * Gets the source file for this read.
     * @param id Unique identifier determining which input file to use.
     * @return The source filename for this read.
     */
    public File getSourceFileForReaderID(final SAMReaderID id) {
        return getReadsDataSource().getSAMFile(id);
    }

    /**
     * Now that all files are open, validate the sequence dictionaries of the reads vs. the reference vrs the reference ordered data (if available).
     *
     * @param reads     Reads data source.
     * @param reference Reference data source.
     * @param rods    a collection of the reference ordered data tracks
     */
    private void validateDataSourcesAgainstReference(SAMDataSource reads, ReferenceSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        if ((reads.isEmpty() && (rods == null || rods.isEmpty())) || reference == null )
            return;

        // Compile a set of sequence names that exist in the reference file.
        SAMSequenceDictionary referenceDictionary = reference.getSequenceDictionary();

        if (!reads.isEmpty()) {
            // Compile a set of sequence names that exist in the BAM files.
            SAMSequenceDictionary readsDictionary = reads.getHeader().getSequenceDictionary();

            if (readsDictionary.isEmpty()) {
                logger.info("Reads file is unmapped. Skipping validation against reference.");
                return;
            }

            // compare the reads to the reference
            SequenceDictionaryUtils.validateDictionaries(logger, getArguments().unsafe, "reads", readsDictionary,
                                                         "reference", referenceDictionary, true, intervals);
        }

        for (ReferenceOrderedDataSource rod : rods)
            IndexDictionaryUtils.validateTrackSequenceDictionary(rod.getName(), rod.getSequenceDictionary(), referenceDictionary, getArguments().unsafe);
    }

    /**
     * Gets a data source for the given set of reads.
     *
     * @param argCollection arguments
     * @param genomeLocParser parser
     * @param refReader reader
     * @return A data source for the given set of reads.
     */
    private SAMDataSource createReadsDataSource(final GATKArgumentCollection argCollection, final GenomeLocParser genomeLocParser,
                                                final IndexedFastaSequenceFile refReader, final Map<String, String> sampleRenameMap) {
        DownsamplingMethod downsamplingMethod = getDownsamplingMethod();

        // Synchronize the method back into the collection so that it shows up when
        // interrogating for the downsampling method during command line recreation.
        setDownsamplingMethod(downsamplingMethod);

        logger.info(downsamplingMethod);

        if (argCollection.removeProgramRecords && argCollection.keepProgramRecords)
            throw new UserException.BadArgumentValue("rpr / kpr", "Cannot enable both options");

        boolean removeProgramRecords = argCollection.removeProgramRecords || walker.getClass().isAnnotationPresent(RemoveProgramRecords.class);

        if (argCollection.keepProgramRecords)
            removeProgramRecords = false;

        final boolean keepReadsInLIBS = walker instanceof ActiveRegionWalker;

        return new SAMDataSource(
                argCollection.referenceFile,
                samReaderIDs,
                threadAllocation,
                argCollection.numberOfBAMFileHandles,
                genomeLocParser,
                argCollection.useOriginalBaseQualities,
                argCollection.strictnessLevel,
                argCollection.readBufferSize,
                downsamplingMethod,
                new ValidationExclusion(Arrays.asList(argCollection.unsafe)),
                filters,
                readTransformers,
                includeReadsWithDeletionAtLoci(),
                argCollection.defaultBaseQualities,
                removeProgramRecords,
                keepReadsInLIBS,
                sampleRenameMap,
                argCollection.intervalArguments.intervalMerging);
    }

    /**
     * Loads a user-provided sample rename map file for use in on-the-fly sample renaming into an in-memory
     * HashMap. This file must consist of lines with two whitespace-separated fields, the second of which
     * may contain whitespace:
     *
     * absolute_path_to_file    new_sample_name
     *
     * The engine will verify that each file contains data from only one sample when the on-the-fly sample
     * renaming feature is being used. Note that this feature works only with bam and vcf files.
     *
     * @param sampleRenameMapFile sample rename map file from which to load data
     * @return a HashMap containing the contents of the map file, with the keys being the input file paths and
     *         the values being the new sample names.
     */
    protected Map<String, String> loadSampleRenameMap( final File sampleRenameMapFile ) {
        logger.info("Renaming samples from input files on-the-fly using mapping file " + sampleRenameMapFile.getAbsolutePath());

        final Map<String, String> sampleRenameMap = new HashMap<>((int)sampleRenameMapFile.length() / 50);

        try {
            for ( final String line : new XReadLines(sampleRenameMapFile) ) {
                final String[] tokens = line.split("\\s+", 2);

                if ( tokens.length != 2 ) {
                    throw new UserException.MalformedFile(sampleRenameMapFile,
                                                          String.format("Encountered a line with %s fields instead of the required 2 fields. Line was: %s",
                                                                        tokens.length, line));
                }

                final File inputFile = new File(tokens[0]);
                final String newSampleName = tokens[1].trim();

                if (newSampleName.contains(VCFConstants.FIELD_SEPARATOR)) {
                    throw new UserException.MalformedFile(sampleRenameMapFile, String.format(
                            "Encountered illegal sample name; sample names may not include the VCF field delimiter (%s).  Sample name: %s; line: %s",
                            VCFConstants.FIELD_SEPARATOR,
                            newSampleName,
                            line
                    ));
                }

                if ( ! inputFile.isAbsolute() ) {
                    throw new UserException.MalformedFile(sampleRenameMapFile, "Input file path not absolute at line: " + line);
                }

                final String inputFilePath = inputFile.getAbsolutePath();

                if ( sampleRenameMap.containsKey(inputFilePath) ) {
                    throw new UserException.MalformedFile(sampleRenameMapFile,
                                                          String.format("Input file %s appears more than once", inputFilePath));
                }

                sampleRenameMap.put(inputFilePath, newSampleName);
            }
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(sampleRenameMapFile, e);
        }

        return sampleRenameMap;
    }


    /**
     * Opens a reference sequence file paired with an index.  Only public for testing purposes
     *
     * @param refFile Handle to a reference sequence file.  Non-null.
     */
    public void setReferenceDataSource(File refFile) {
        this.referenceDataSource = new ReferenceDataSource(refFile);
        genomeLocParser = new GenomeLocParser(referenceDataSource.getReference());
    }

    /**
     * Open the reference-ordered data sources.
     *
     * @param referenceMetaDataFiles collection of RMD descriptors to load and validate.
     * @param sequenceDictionary GATK-wide sequnce dictionary to use for validation.
     * @param genomeLocParser to use when creating and validating GenomeLocs.
     * @param validationExclusionType potentially indicate which validations to include / exclude.
     * @param sampleRenameMap map of file -> new sample name used when doing on-the-fly sample renaming
     *
     * @return A list of reference-ordered data sources.
     */
    private List<ReferenceOrderedDataSource> getReferenceOrderedDataSources(final Collection<RMDTriplet> referenceMetaDataFiles,
                                                                            final SAMSequenceDictionary sequenceDictionary,
                                                                            final GenomeLocParser genomeLocParser,
                                                                            final ValidationExclusion.TYPE validationExclusionType,
                                                                            final Map<String, String> sampleRenameMap) {
        final RMDTrackBuilder builder = new RMDTrackBuilder(sequenceDictionary,genomeLocParser, validationExclusionType,
                                                            getArguments().disableAutoIndexCreationAndLockingWhenReadingRods,
                                                            sampleRenameMap);

        final List<ReferenceOrderedDataSource> dataSources = new ArrayList<ReferenceOrderedDataSource>();
        for (RMDTriplet fileDescriptor : referenceMetaDataFiles)
            dataSources.add(new ReferenceOrderedDataSource(fileDescriptor,
                                                           builder,
                                                           sequenceDictionary,
                                                           genomeLocParser,
                                                           flashbackData()));

        return dataSources;
    }

    /**
     * Returns the SAM File Header from the input reads' data source file
     * @return the SAM File Header from the input reads' data source file
     */
    public SAMFileHeader getSAMFileHeader() {
        return readsDataSource.getHeader();
    }

    public boolean lenientVCFProcessing() {
        return ValidationExclusion.lenientVCFProcessing(argCollection.unsafe);
    }

    /**
     * Returns the unmerged SAM file header for an individual reader.
     * @param reader The reader.
     * @return Header for that reader or null if not available.
     */
    public SAMFileHeader getSAMFileHeader(SAMReaderID reader) {
        return readsDataSource == null ? null : readsDataSource.getHeader(reader);
    }

    /**
     * Returns an ordered list of the unmerged SAM file headers known to this engine.
     * @return list of header for each input SAM file, in command line order
     */
    public List<SAMFileHeader> getSAMFileHeaders() {
        final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
        for ( final SAMReaderID id : getReadsDataSource().getReaderIDs() ) {
            headers.add(getReadsDataSource().getHeader(id));
        }
        return headers;
    }

    /**
     * Gets the master sequence dictionary for this GATK engine instance
     * @return a never-null dictionary listing all of the contigs known to this engine instance
     */
    public SAMSequenceDictionary getMasterSequenceDictionary() {
        return getReferenceDataSource().getReference().getSequenceDictionary();
    }

    /**
     * Returns data source object encapsulating all essential info and handlers used to traverse
     * reads; header merger, individual file readers etc can be accessed through the returned data source object.
     *
     * @return the reads data source
     */
    public SAMDataSource getReadsDataSource() {
        return this.readsDataSource;
    }

    /**
     * Sets the collection of GATK main application arguments.
     *
     * @param argCollection the GATK argument collection
     */
    public void setArguments(GATKArgumentCollection argCollection) {
        this.argCollection = argCollection;
    }

    /**
     * Gets the collection of GATK main application arguments.
     *
     * @return the GATK argument collection
     */
    public GATKArgumentCollection getArguments() {
        return this.argCollection;
    }

    /**
     * Get the list of intervals passed to the engine.
     * @return List of intervals, or null if no intervals are in use
     */
    public GenomeLocSortedSet getIntervals() {
        return this.intervals;
    }

    /**
     * Get the list of regions of the genome being processed.  If the user
     * requested specific intervals, return those, otherwise return regions
     * corresponding to the entire genome.  Never returns null.
     *
     * @return a non-null set of intervals being processed
     */
    @Ensures("result != null")
    public GenomeLocSortedSet getRegionsOfGenomeBeingProcessed() {
        if ( getIntervals() == null )
            // if we don't have any intervals defined, create intervals from the reference itself
            return GenomeLocSortedSet.createSetFromSequenceDictionary(getReferenceDataSource().getReference().getSequenceDictionary());
        else
            return getIntervals();
    }

    /**
     * Gets the list of filters employed by this engine.
     * @return Collection of filters (actual instances) used by this engine.
     */
    public Collection<ReadFilter> getFilters() {
        return this.filters;
    }

    /**
     * Sets the list of filters employed by this engine.
     * @param filters Collection of filters (actual instances) used by this engine.
     */
    public void setFilters(Collection<ReadFilter> filters) {
        this.filters = filters;
    }

    /**
     * Gets the filter manager for this engine.
     * @return filter manager for this engine.
     */
    protected FilterManager getFilterManager() {
        return filterManager;
    }

    /**
     * Gets the input sources for this engine.
     * @return input sources for this engine.
     */
    protected Map<ArgumentSource, Object> getInputs() {
        return inputs;
    }

    /**
     * Gets the output stubs for this engine.
     * @return output stubs for this engine.
     */
    protected Collection<Stub<?>> getOutputs() {
        return outputs;
    }

    /**
     * Returns data source objects encapsulating all rod data;
     * individual rods can be accessed through the returned data source objects.
     *
     * @return the rods data sources, never {@code null}.
     */
    public List<ReferenceOrderedDataSource> getRodDataSources() {
        return this.rodDataSources;
    }

    /**
     * Gets cumulative metrics about the entire run to this point.
     * Returns a clone of this snapshot in time.
     * @return cumulative metrics about the entire run at this point.  ReadMetrics object is a unique instance and is
     *         owned by the caller; the caller can do with the object what they wish.
     */
    public ReadMetrics getCumulativeMetrics() {
        // todo -- probably shouldn't be lazy
        if ( cumulativeMetrics == null )
            cumulativeMetrics = readsDataSource == null ? new ReadMetrics() : readsDataSource.getCumulativeReadMetrics();
        return cumulativeMetrics;
    }

    /**
     * Return the global ThreadEfficiencyMonitor, if there is one
     *
     * @return the monitor, or null if none is active
     */
    public ThreadEfficiencyMonitor getThreadEfficiencyMonitor() {
        return threadEfficiencyMonitor;
    }

    // -------------------------------------------------------------------------------------
    //
    // code for working with Samples database
    //
    // -------------------------------------------------------------------------------------

    public SampleDB getSampleDB() {
        return this.sampleDB;
    }

    public Map<String,String> getApproximateCommandLineArguments(Object... argumentProviders) {
        return CommandLineUtils.getApproximateCommandLineArguments(parsingEngine,argumentProviders);
    }

    public String createApproximateCommandLineArgumentString(Object... argumentProviders) {
        return CommandLineUtils.createApproximateCommandLineArgumentString(parsingEngine,argumentProviders);
    }

    // -------------------------------------------------------------------------------------
    //
    // code for working with progress meter
    //
    // -------------------------------------------------------------------------------------

    /**
     * Register the global progress meter with this engine
     *
     * Calling this function more than once will result in an IllegalStateException
     *
     * @param meter a non-null progress meter
     */
    public void registerProgressMeter(final ProgressMeter meter) {
        if ( meter == null ) throw new IllegalArgumentException("Meter cannot be null");
        if ( progressMeter != null ) throw new IllegalStateException("Progress meter already set");

        progressMeter = meter;
    }

    /**
     * Get the progress meter being used by this engine.  May be null if no meter has been registered yet
     * @return a potentially null pointer to the progress meter
     */
    public ProgressMeter getProgressMeter() {
        return progressMeter;
    }

    /**
     * Does the current runtime in unit exceed the runtime limit, if one has been provided?
     *
     * @return false if not limit was requested or if runtime <= the limit, true otherwise
     */
    public boolean exceedsRuntimeLimit() {
        if ( progressMeter == null )
            // not yet initialized or not set because of testing
            return false;

        if ( getArguments().maxRuntime == NO_RUNTIME_LIMIT )
            return false;
        else {  
            final long runtime = progressMeter.getRuntimeInNanosecondsUpdatedPeriodically();
            if ( runtime < 0 ) throw new IllegalArgumentException("runtime must be >= 0 but got " + runtime);
            final long maxRuntimeNano = getRuntimeLimitInNanoseconds();
            return runtime > maxRuntimeNano;
        }
    }

    /**
     * @return the runtime limit in nanoseconds, or -1 if no limit was specified
     */
    public long getRuntimeLimitInNanoseconds() {
        return runtimeLimitInNanoseconds;
    }

    /**
     * Setup the runtime limits for this engine, updating the runtimeLimitInNanoseconds
     * as appropriate
     *
     * @param args the GATKArgumentCollection to retrieve our runtime limits from
     */
    private void setupRuntimeLimits(final GATKArgumentCollection args) {
        if ( args.maxRuntime == NO_RUNTIME_LIMIT )
            runtimeLimitInNanoseconds = -1;
        else if (args.maxRuntime < 0 )
            throw new UserException.BadArgumentValue("maxRuntime", "must be >= 0 or == -1 (meaning no limit) but received negative value " + args.maxRuntime);
        else {
            runtimeLimitInNanoseconds = TimeUnit.NANOSECONDS.convert(args.maxRuntime, args.maxRuntimeUnits);
        }
    }

    /**
     * Returns the sample list including all samples.
     * @return never {@code null}.
     */
    public SampleList getSampleList() {
        return new IndexedSampleList(getSampleDB().getSampleNames());
    }

    /**
     * Returns the sample list including samples in read inputs.
     * @return never {@code null}.
     */
    public SampleList getReadSampleList() {
        return new IndexedSampleList(ReadUtils.getSAMFileSamples(getSAMFileHeader()));
    }
}
