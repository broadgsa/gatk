/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.reads.*;
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.samples.SampleDB;
import org.broadinstitute.sting.gatk.executive.MicroScheduler;
import org.broadinstitute.sting.gatk.filters.FilterManager;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.filters.ReadGroupBlackListFilter;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.io.stubs.Stub;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.samples.SampleDBBuilder;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;

import java.io.File;
import java.util.*;

/**
 * A GenomeAnalysisEngine that runs a specified walker.
 */
public class GenomeAnalysisEngine {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(GenomeAnalysisEngine.class);

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
    private SampleDB sampleDB = null;

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
     * A currently hacky unique name for this GATK instance
     */
    private String myName = "GATK_" + Math.abs(getRandomGenerator().nextInt());

    /**
     * our walker manager
     */
    private final WalkerManager walkerManager = new WalkerManager();

    private Walker<?, ?> walker;

    public void setWalker(Walker<?, ?> walker) {
        this.walker = walker;
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
     * Set the reference metadata files to use for this traversal.
     * @param referenceMetaDataFiles Collection of files and descriptors over which to traverse.
     */
    public void setReferenceMetaDataFiles(Collection<RMDTriplet> referenceMetaDataFiles) {
        this.referenceMetaDataFiles = referenceMetaDataFiles;
    }

    /**
     *  Static random number generator and seed.
     */
    private static final long GATK_RANDOM_SEED = 47382911L;
    private static Random randomGenerator = new Random(GATK_RANDOM_SEED);

    public static Random getRandomGenerator() { return randomGenerator; }
    public static void resetRandomGenerator() { randomGenerator.setSeed(GATK_RANDOM_SEED); }
    public static void resetRandomGenerator(long seed) { randomGenerator.setSeed(seed); }
    /**
     * Actually run the GATK with the specified walker.
     *
     * @return the value of this traversal.
     */
    public Object execute() {
        //HeapSizeMonitor monitor = new HeapSizeMonitor();
        //monitor.start();
        setStartTime(new java.util.Date());

        // validate our parameters
        if (this.getArguments() == null) {
            throw new ReviewedStingException("The GATKArgumentCollection passed to GenomeAnalysisEngine can not be null.");
        }

        // validate our parameters              
        if (this.walker == null)
            throw new ReviewedStingException("The walker passed to GenomeAnalysisEngine can not be null.");

        if (this.getArguments().nonDeterministicRandomSeed)
            resetRandomGenerator(System.currentTimeMillis());

        // Prepare the data for traversal.
        initializeDataSources();

        // initialize sampleDB
        initializeSampleDB();

        // initialize and validate the interval list
        initializeIntervals();
        validateSuppliedIntervals();

        // our microscheduler, which is in charge of running everything
        MicroScheduler microScheduler = createMicroscheduler();

        // create temp directories as necessary
        initializeTempDirectory();

        // create the output streams                     "
        initializeOutputStreams(microScheduler.getOutputTracker());

        ShardStrategy shardStrategy = getShardStrategy(readsDataSource,microScheduler.getReference(),intervals);

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
        return walkerManager.createByName(walkerName);
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
        final List<ReadFilter> filters = WalkerManager.getReadFilters(walker,this.getFilterManager());
        if (this.getArguments().readGroupBlackList != null && this.getArguments().readGroupBlackList.size() > 0)
            filters.add(new ReadGroupBlackListFilter(this.getArguments().readGroupBlackList));
        for(final String filterName: this.getArguments().readFilters)
            filters.add(this.getFilterManager().createByName(filterName));
        return Collections.unmodifiableList(filters);
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

        return MicroScheduler.create(this,walker,this.getReadsDataSource(),this.getReferenceDataSource().getReference(),this.getRodDataSources(),this.getArguments().numberOfThreads);
    }

    protected DownsamplingMethod getDownsamplingMethod() {
        GATKArgumentCollection argCollection = this.getArguments();
        DownsamplingMethod method;
        if(argCollection.getDownsamplingMethod() != null)
            method = argCollection.getDownsamplingMethod();
        else if(WalkerManager.getDownsamplingMethod(walker) != null)
            method = WalkerManager.getDownsamplingMethod(walker);
        else
            method = GATKArgumentCollection.getDefaultDownsamplingMethod();
        return method;
    }

    protected void setDownsamplingMethod(DownsamplingMethod method) {
        argCollection.setDownsamplingMethod(method);
    }

    public BAQ.QualityMode getWalkerBAQQualityMode()         { return WalkerManager.getBAQQualityMode(walker); }
    public BAQ.ApplicationTime getWalkerBAQApplicationTime() { return WalkerManager.getBAQApplicationTime(walker); }    

    protected boolean generateExtendedEvents() {
        return walker.generateExtendedEvents();
    }

    protected boolean includeReadsWithDeletionAtLoci() {
        return walker.includeReadsWithDeletionAtLoci();
    }

    /**
     * Verifies that the supplied set of reads files mesh with what the walker says it requires,
     * and also makes sure that there were no duplicate SAM files specified on the command line.
     */
    protected void validateSuppliedReads() {
        GATKArgumentCollection arguments = this.getArguments();
        // Check what the walker says is required against what was provided on the command line.
        if (WalkerManager.isRequired(walker, DataSource.READS) && (arguments.samFiles == null || arguments.samFiles.size() == 0))
            throw new ArgumentException("Walker requires reads but none were provided.");

        // Check what the walker says is allowed against what was provided on the command line.
        if ((arguments.samFiles != null && arguments.samFiles.size() > 0) && !WalkerManager.isAllowed(walker, DataSource.READS))
            throw new ArgumentException("Walker does not allow reads but reads were provided.");

        // Make sure no SAM files were specified multiple times by the user.
        checkForDuplicateSamFiles();
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
            throw new ArgumentException("The following BAM files appear multiple times in the list of input files: " +
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
    }

    /**
     * Get the sharding strategy given a driving data source.
     *
     * @param readsDataSource readsDataSource
     * @param drivingDataSource Data on which to shard.
     * @param intervals intervals
     * @return the sharding strategy
     */
    protected ShardStrategy getShardStrategy(SAMDataSource readsDataSource, ReferenceSequenceFile drivingDataSource, GenomeLocSortedSet intervals) {
        ValidationExclusion exclusions = (readsDataSource != null ? readsDataSource.getReadsInfo().getValidationExclusionList() : null);
        ReferenceDataSource referenceDataSource = this.getReferenceDataSource();
        // Use monolithic sharding if no index is present.  Monolithic sharding is always required for the original
        // sharding system; it's required with the new sharding system only for locus walkers.
        if(readsDataSource != null && !readsDataSource.hasIndex() ) { 
            if(!exclusions.contains(ValidationExclusion.TYPE.ALLOW_UNINDEXED_BAM))
                throw new UserException.CommandLineException("Cannot process the provided BAM file(s) because they were not indexed.  The GATK does offer limited processing of unindexed BAMs in --unsafe mode, but this GATK feature is currently unsupported.");
            if(intervals != null && !argCollection.allowIntervalsWithUnindexedBAM)
                throw new UserException.CommandLineException("Cannot perform interval processing when reads are present but no index is available.");

            Shard.ShardType shardType;
            if(walker instanceof LocusWalker) {
                if (readsDataSource.getSortOrder() != SAMFileHeader.SortOrder.coordinate)
                    throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.coordinate, "Locus walkers can only traverse coordinate-sorted data.  Please resort your input BAM file(s) or set the Sort Order tag in the header appropriately.");
                shardType = Shard.ShardType.LOCUS;
            }
            else if(walker instanceof ReadWalker || walker instanceof DuplicateWalker || walker instanceof ReadPairWalker)
                shardType = Shard.ShardType.READ;
            else
                throw new UserException.CommandLineException("The GATK cannot currently process unindexed BAM files");

            List<GenomeLoc> region;
            if(intervals != null)
                region = intervals.toList();
            else {
                region = new ArrayList<GenomeLoc>();
                for(SAMSequenceRecord sequenceRecord: drivingDataSource.getSequenceDictionary().getSequences())
                    region.add(getGenomeLocParser().createGenomeLoc(sequenceRecord.getSequenceName(),1,sequenceRecord.getSequenceLength()));
            }

            return new MonolithicShardStrategy(getGenomeLocParser(), readsDataSource,shardType,region);
        }

        ShardStrategy shardStrategy;
        ShardStrategyFactory.SHATTER_STRATEGY shardType;

        long SHARD_SIZE = 100000L;

        if (walker instanceof LocusWalker) {
            if (walker instanceof RodWalker) SHARD_SIZE *= 1000;

            if (intervals != null && !intervals.isEmpty()) {
                if (readsDataSource == null)
                    throw new IllegalArgumentException("readsDataSource is null");
                if(!readsDataSource.isEmpty() && readsDataSource.getSortOrder() != SAMFileHeader.SortOrder.coordinate)
                    throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.coordinate, "Locus walkers can only traverse coordinate-sorted data.  Please resort your input BAM file(s) or set the Sort Order tag in the header appropriately.");

                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        referenceDataSource.getReference(),
                        ShardStrategyFactory.SHATTER_STRATEGY.LOCUS_EXPERIMENTAL,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE,
                        getGenomeLocParser(),
                        intervals);
            } else
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        referenceDataSource.getReference(),
                        ShardStrategyFactory.SHATTER_STRATEGY.LOCUS_EXPERIMENTAL,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE,getGenomeLocParser());
        } else if (walker instanceof ReadWalker ||
                walker instanceof DuplicateWalker) {
            shardType = ShardStrategyFactory.SHATTER_STRATEGY.READS_EXPERIMENTAL;

            if (intervals != null && !intervals.isEmpty()) {
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        referenceDataSource.getReference(),
                        shardType,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE,
                        getGenomeLocParser(),
                        intervals);
            } else {
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        referenceDataSource.getReference(),
                        shardType,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE,
                        getGenomeLocParser());
            }
        } else if (walker instanceof ReadPairWalker) {
            if(readsDataSource != null && readsDataSource.getSortOrder() != SAMFileHeader.SortOrder.queryname)
                throw new UserException.MissortedBAM(SAMFileHeader.SortOrder.queryname, "Read pair walkers can only walk over query name-sorted data.  Please resort your input BAM file.");
            if(intervals != null && !intervals.isEmpty())
                throw new UserException.CommandLineException("Pairs traversal cannot be used in conjunction with intervals.");

            shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                    referenceDataSource.getReference(),
                    ShardStrategyFactory.SHATTER_STRATEGY.READS_EXPERIMENTAL,
                    drivingDataSource.getSequenceDictionary(),
                    SHARD_SIZE,
                    getGenomeLocParser());
        } else
            throw new ReviewedStingException("Unable to support walker of type" + walker.getClass().getName());

        return shardStrategy;
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
    private void initializeOutputStreams(OutputTracker outputTracker) {
        for (Map.Entry<ArgumentSource, Object> input : getInputs().entrySet())
            outputTracker.addInput(input.getKey(), input.getValue());
        for (Stub<?> stub : getOutputs())
            outputTracker.addOutput(stub);

        outputTracker.prepareWalker(walker, getArguments().strictnessLevel);
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

        // return if no interval arguments at all
        if ( argCollection.intervals == null && argCollection.excludeIntervals == null )
            return;

        // Note that the use of '-L all' is no longer supported.

        // if include argument isn't given, create new set of all possible intervals
        GenomeLocSortedSet includeSortedSet = (argCollection.intervals == null ?
            GenomeLocSortedSet.createSetFromSequenceDictionary(this.referenceDataSource.getReference().getSequenceDictionary()) :
            loadIntervals(argCollection.intervals, argCollection.intervalSetRule));

        // if no exclude arguments, can return parseIntervalArguments directly
        if ( argCollection.excludeIntervals == null )
            intervals = includeSortedSet;

        // otherwise there are exclude arguments => must merge include and exclude GenomeLocSortedSets
        else {
            GenomeLocSortedSet excludeSortedSet = loadIntervals(argCollection.excludeIntervals, IntervalSetRule.UNION);
            intervals = includeSortedSet.subtractRegions(excludeSortedSet);

            // logging messages only printed when exclude (-XL) arguments are given
            long toPruneSize = includeSortedSet.coveredSize();
            long toExcludeSize = excludeSortedSet.coveredSize();
            long intervalSize = intervals.coveredSize();
            logger.info(String.format("Initial include intervals span %d loci; exclude intervals span %d loci", toPruneSize, toExcludeSize));
            logger.info(String.format("Excluding %d loci from original intervals (%.2f%% reduction)",
                    toPruneSize - intervalSize, (toPruneSize - intervalSize) / (0.01 * toPruneSize)));
        }
    }

    /**
     * Loads the intervals relevant to the current execution
     * @param argList  argument bindings; might include filenames, intervals in samtools notation, or a combination of the above
     * @param rule     interval merging rule
     * @return A sorted, merged list of all intervals specified in this arg list.
     */
    protected GenomeLocSortedSet loadIntervals( List<IntervalBinding<Feature>> argList, IntervalSetRule rule ) {

        List<GenomeLoc> allIntervals = new ArrayList<GenomeLoc>(0);
        for ( IntervalBinding intervalBinding : argList ) {
            List<GenomeLoc> intervals = intervalBinding.getIntervals(this);

            if ( intervals.isEmpty() ) {
                logger.warn("The interval file " + intervalBinding.getSource() + " contains no intervals that could be parsed.");
            }

            allIntervals = IntervalUtils.mergeListsBySetOperator(intervals, allIntervals, rule);
        }

        return IntervalUtils.sortAndMergeIntervals(genomeLocParser, allIntervals, argCollection.intervalMerging);
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

        // TODO -- REMOVE ME
        BAQ.DEFAULT_GOP = argCollection.BAQGOP;

        validateSuppliedReference();
        setReferenceDataSource(argCollection.referenceFile);

        validateSuppliedReads();
        readsDataSource = createReadsDataSource(argCollection,genomeLocParser,referenceDataSource.getReference());

        for (ReadFilter filter : filters)
            filter.initialize(this);

        // set the sequence dictionary of all of Tribble tracks to the sequence dictionary of our reference
        rodDataSources = getReferenceOrderedDataSources(referenceMetaDataFiles,referenceDataSource.getReference().getSequenceDictionary(),genomeLocParser,argCollection.unsafe);
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
     * @param manager manager
     */
    private void validateSourcesAgainstReference(SAMDataSource reads, ReferenceSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, RMDTrackBuilder manager) {
        if ((reads.isEmpty() && (rods == null || rods.isEmpty())) || reference == null )
            return;

        // Compile a set of sequence names that exist in the reference file.
        SAMSequenceDictionary referenceDictionary = reference.getSequenceDictionary();

        if (!reads.isEmpty()) {
            // Compile a set of sequence names that exist in the BAM files.
            SAMSequenceDictionary readsDictionary = reads.getHeader().getSequenceDictionary();

            if (readsDictionary.size() == 0) {
                logger.info("Reads file is unmapped.  Skipping validation against reference.");
                return;
            }

            // compare the reads to the reference
            SequenceDictionaryUtils.validateDictionaries(logger, getArguments().unsafe, "reads", readsDictionary, "reference", referenceDictionary);
        }

        for (ReferenceOrderedDataSource rod : rods)
            manager.validateTrackSequenceDictionary(rod.getName(),rod.getSequenceDictionary(),referenceDictionary);
    }

    /**
     * Gets a data source for the given set of reads.
     *
     * @param argCollection arguments
     * @param genomeLocParser parser
     * @param refReader reader
     * @return A data source for the given set of reads.
     */
    private SAMDataSource createReadsDataSource(GATKArgumentCollection argCollection, GenomeLocParser genomeLocParser, IndexedFastaSequenceFile refReader) {
        DownsamplingMethod method = getDownsamplingMethod();

        // Synchronize the method back into the collection so that it shows up when
        // interrogating for the downsample method during command line recreation.
        setDownsamplingMethod(method);

        if ( getWalkerBAQApplicationTime() == BAQ.ApplicationTime.FORBIDDEN && argCollection.BAQMode != BAQ.CalculationMode.OFF)
            throw new UserException.BadArgumentValue("baq", "Walker cannot accept BAQ'd base qualities, and yet BAQ mode " + argCollection.BAQMode + " was requested.");

        return new SAMDataSource(
                samReaderIDs,
                genomeLocParser,
                argCollection.useOriginalBaseQualities,
                argCollection.strictnessLevel,
                argCollection.readBufferSize,
                method,
                new ValidationExclusion(Arrays.asList(argCollection.unsafe)),
                filters,
                includeReadsWithDeletionAtLoci(),
                generateExtendedEvents(),
                getWalkerBAQApplicationTime() == BAQ.ApplicationTime.ON_INPUT ? argCollection.BAQMode : BAQ.CalculationMode.OFF,
                getWalkerBAQQualityMode(),
                refReader,
                argCollection.defaultBaseQualities,
                !argCollection.disableLowMemorySharding);
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
     *
     * @return A list of reference-ordered data sources.
     */
    private List<ReferenceOrderedDataSource> getReferenceOrderedDataSources(Collection<RMDTriplet> referenceMetaDataFiles,
                                                                            SAMSequenceDictionary sequenceDictionary,
                                                                            GenomeLocParser genomeLocParser,
                                                                            ValidationExclusion.TYPE validationExclusionType) {
        RMDTrackBuilder builder = new RMDTrackBuilder(sequenceDictionary,genomeLocParser,validationExclusionType);

        List<ReferenceOrderedDataSource> dataSources = new ArrayList<ReferenceOrderedDataSource>();
        for (RMDTriplet fileDescriptor : referenceMetaDataFiles)
            dataSources.add(new ReferenceOrderedDataSource(fileDescriptor,
                                                           builder,
                                                           sequenceDictionary,
                                                           genomeLocParser,
                                                           flashbackData()));

        // validation: check to make sure everything the walker needs is present, and that all sequence dictionaries match.
        validateSourcesAgainstReference(readsDataSource, referenceDataSource.getReference(), dataSources, builder);

        return dataSources;
    }

    /**
     * Returns the SAM File Header from the input reads' data source file
     * @return the SAM File Header from the input reads' data source file
     */
    public SAMFileHeader getSAMFileHeader() {
        return readsDataSource.getHeader();
    }

    /**
     * Returns the unmerged SAM file header for an individual reader.
     * @param reader The reader.
     * @return Header for that reader.
     */
    public SAMFileHeader getSAMFileHeader(SAMReaderID reader) {
        return readsDataSource.getHeader(reader);
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
     * @return the rods data sources
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
        return readsDataSource == null ? null : readsDataSource.getCumulativeReadMetrics();
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
    

}
