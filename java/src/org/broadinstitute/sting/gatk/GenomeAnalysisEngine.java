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

import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.ArgumentException;
import org.broadinstitute.sting.commandline.ArgumentSource;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.commandline.ParsingEngine;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.datasources.sample.SampleDataSource;
import org.broadinstitute.sting.gatk.datasources.shards.MonolithicShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.gatk.executive.MicroScheduler;
import org.broadinstitute.sting.gatk.filters.FilterManager;
import org.broadinstitute.sting.gatk.filters.ReadGroupBlackListFilter;
import org.broadinstitute.sting.gatk.filters.SamRecordHeaderFilter;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.io.stubs.Stub;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.RMDIntervalGenerator;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.SequenceDictionaryUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
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
    private SampleDataSource sampleDataSource = null;

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
    private Collection<SamRecordFilter> filters;

    /**
     * A currently hacky unique name for this GATK instance
     */
    private String myName = "GATK_" + Math.abs(new Random().nextInt());

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

        // Prepare the data for traversal.
        initializeDataSources();

        // initialize and validate the interval list
        initializeIntervals();
        validateSuppliedIntervals();

        // our microscheduler, which is in charge of running everything
        MicroScheduler microScheduler = createMicroscheduler();

        // create temp directories as necessary
        initializeTempDirectory();

        // create the output streams                     "
        initializeOutputStreams(microScheduler.getOutputTracker());

        ShardStrategy shardStrategy = getShardStrategy(microScheduler.getReference());

        // execute the microscheduler, storing the results
        Object result =  microScheduler.execute(this.walker, shardStrategy);

        //monitor.stop();
        //logger.info(String.format("Maximum heap size consumed: %d",monitor.getMaxMemoryUsed()));

        return result;
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
    public Collection<SamRecordFilter> createFilters() {
        Set<SamRecordFilter> filters = new HashSet<SamRecordFilter>();
        filters.addAll(WalkerManager.getReadFilters(walker,this.getFilterManager()));
        if (this.getArguments().readGroupBlackList != null && this.getArguments().readGroupBlackList.size() > 0)
            filters.add(new ReadGroupBlackListFilter(this.getArguments().readGroupBlackList));
        for(String filterName: this.getArguments().readFilters)
            filters.add(this.getFilterManager().createByName(filterName));
        return Collections.unmodifiableSet(filters);
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
            method = argCollection.getDefaultDownsamplingMethod();
        return method;
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
     * Verifies that the supplied set of reads files mesh with what the walker says it requires.
     */
    protected void validateSuppliedReads() {
        GATKArgumentCollection arguments = this.getArguments();
        // Check what the walker says is required against what was provided on the command line.
        if (WalkerManager.isRequired(walker, DataSource.READS) && (arguments.samFiles == null || arguments.samFiles.size() == 0))
            throw new ArgumentException("Walker requires reads but none were provided.");

        // Check what the walker says is allowed against what was provided on the command line.
        if ((arguments.samFiles != null && arguments.samFiles.size() > 0) && !WalkerManager.isAllowed(walker, DataSource.READS))
            throw new ArgumentException("Walker does not allow reads but reads were provided.");
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

    /**
     * Verifies that all required reference-ordered data has been supplied, and any reference-ordered data that was not
     * 'allowed' is still present.
     *
     * @param rods   Reference-ordered data to load.
     */
    protected void validateSuppliedReferenceOrderedData(List<ReferenceOrderedDataSource> rods) {
        // Check to make sure that all required metadata is present.
        List<RMD> allRequired = WalkerManager.getRequiredMetaData(walker);
        for (RMD required : allRequired) {
            boolean found = false;
            for (ReferenceOrderedDataSource rod : rods) {
                if (rod.matchesNameAndRecordType(required.name(), required.type()))
                    found = true;
            }
            if (!found)
                throw new ArgumentException(String.format("Walker requires reference metadata to be supplied named '%s' of type '%s', but this metadata was not provided.  " +
                                                          "Please supply the specified metadata file.", required.name(), required.type().getSimpleName()));
        }

        // Check to see that no forbidden rods are present.
        for (ReferenceOrderedDataSource rod : rods) {
            if (!WalkerManager.isAllowed(walker, rod))
                throw new ArgumentException(String.format("Walker of type %s does not allow access to metadata: %s", walker.getClass(), rod.getName()));
        }
    }

    protected void validateSuppliedIntervals() {
        // Only read walkers support '-L unmapped' intervals.  Trap and validate any other instances of -L unmapped.
        if(!(walker instanceof ReadWalker)) {
            GenomeLocSortedSet intervals = getIntervals();
            if(intervals != null && getIntervals().contains(GenomeLoc.UNMAPPED))
                throw new ArgumentException("Interval list specifies unmapped region.  Only read walkers may include the unmapped region.");
        }
    }

    /**
     * Get the sharding strategy given a driving data source.
     *
     * @param drivingDataSource Data on which to shard.
     * @return the sharding strategy
     */
    protected ShardStrategy getShardStrategy(ReferenceSequenceFile drivingDataSource) {
        GenomeLocSortedSet intervals = this.getIntervals();
        SAMDataSource readsDataSource = this.getReadsDataSource();
        ValidationExclusion exclusions = (readsDataSource != null ? readsDataSource.getReadsInfo().getValidationExclusionList() : null);
        ReferenceDataSource referenceDataSource = this.getReferenceDataSource();
        // Use monolithic sharding if no index is present.  Monolithic sharding is always required for the original
        // sharding system; it's required with the new sharding system only for locus walkers.
        if(readsDataSource != null && !readsDataSource.hasIndex() ) { 
            if(!exclusions.contains(ValidationExclusion.TYPE.ALLOW_UNINDEXED_BAM))
                throw new UserException.CommandLineException("Cannot process the provided BAM file(s) because they were not indexed.  The GATK does offer limited processing of unindexed BAMs in --unsafe mode, but this GATK feature is currently unsupported.");
            if(intervals != null)
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

            return new MonolithicShardStrategy(readsDataSource,shardType,region);
        }

        ShardStrategy shardStrategy = null;
        ShardStrategyFactory.SHATTER_STRATEGY shardType;

        long SHARD_SIZE = 100000L;

        if (walker instanceof LocusWalker) {
            if (walker instanceof RodWalker) SHARD_SIZE *= 1000;

            if (intervals != null && !intervals.isEmpty()) {
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
        tempDir.mkdirs();
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
        if ((argCollection.intervals == null) && (argCollection.excludeIntervals == null) && (argCollection.RODToInterval == null))
            return;

        // if '-L all' was specified, verify that it was the only -L specified and return if so.
        if(argCollection.intervals != null) {
            for(String interval: argCollection.intervals) {
                if(interval.trim().equals("all")) {
                    if(argCollection.intervals.size() > 1)
                        throw new UserException("'-L all' was specified along with other intervals or interval lists; the GATK cannot combine '-L all' with other intervals.");

                    // '-L all' was specified and seems valid.  Return.
                    return;
                }
            }
        }

        // if include argument isn't given, create new set of all possible intervals
        GenomeLocSortedSet includeSortedSet = (argCollection.intervals == null && argCollection.RODToInterval == null ?
                GenomeLocSortedSet.createSetFromSequenceDictionary(this.referenceDataSource.getReference().getSequenceDictionary()) :
                loadIntervals(argCollection.intervals,
                        argCollection.intervalMerging,
                        genomeLocParser.mergeIntervalLocations(checkRODToIntervalArgument(),argCollection.intervalMerging)));

        // if no exclude arguments, can return parseIntervalArguments directly
        if (argCollection.excludeIntervals == null)
            intervals = includeSortedSet;

            // otherwise there are exclude arguments => must merge include and exclude GenomeLocSortedSets
        else {
            GenomeLocSortedSet excludeSortedSet = loadIntervals(argCollection.excludeIntervals, argCollection.intervalMerging, null);
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
     * @param argList String representation of arguments; might include 'all', filenames, intervals in samtools
     *                notation, or a combination of the
     * @param mergingRule Technique to use when merging interval data.
     * @param additionalIntervals a list of additional intervals to add to the returned set.  Can be null.
     * @return A sorted, merged list of all intervals specified in this arg list.
     */
    private GenomeLocSortedSet loadIntervals(List<String> argList,
                                             IntervalMergingRule mergingRule,
                                             List<GenomeLoc> additionalIntervals) {

        return IntervalUtils.sortAndMergeIntervals(genomeLocParser,IntervalUtils.mergeListsBySetOperator(additionalIntervals,
                IntervalUtils.parseIntervalArguments(genomeLocParser,argList,
                        this.getArguments().unsafe != ValidationExclusion.TYPE.ALLOW_EMPTY_INTERVAL_LIST),
                argCollection.BTIMergeRule),
                mergingRule);
    }

    /**
     * if we have a ROD specified as a 'rodToIntervalTrackName', convert its records to RODs
     * @return ROD intervals as GenomeLocs
     */
    private List<GenomeLoc> checkRODToIntervalArgument() {
        Map<String, ReferenceOrderedDataSource> rodNames = RMDIntervalGenerator.getRMDTrackNames(rodDataSources);
        // Do we have any RODs that overloaded as interval lists with the 'rodToIntervalTrackName' flag?
        List<GenomeLoc> ret = new ArrayList<GenomeLoc>();
        if (rodNames != null && argCollection.RODToInterval != null) {
            String rodName = argCollection.RODToInterval;

            // check to make sure we have a rod of that name
            if (!rodNames.containsKey(rodName))
                throw new UserException.CommandLineException("--rodToIntervalTrackName (-BTI) was passed the name '"+rodName+"', which wasn't given as a ROD name in the -B option");

            for (String str : rodNames.keySet())
                if (str.equals(rodName)) {
                    logger.info("Adding interval list from track (ROD) named " + rodName);
                    RMDIntervalGenerator intervalGenerator = new RMDIntervalGenerator(rodNames.get(str));
                    ret.addAll(intervalGenerator.toGenomeLocList());
                }
        }
        return ret;
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

    protected void initializeDataSources() {
        logger.info("Strictness is " + argCollection.strictnessLevel);

        // TODO -- REMOVE ME
        BAQ.DEFAULT_GOP = argCollection.BAQGOP;

        validateSuppliedReference();
        setReferenceDataSource(argCollection.referenceFile);

        validateSuppliedReads();
        readsDataSource = createReadsDataSource(genomeLocParser, referenceDataSource.getReference());

        sampleDataSource = new SampleDataSource(getSAMFileHeader(), argCollection.sampleFiles);

        for (SamRecordFilter filter : filters)
            if (filter instanceof SamRecordHeaderFilter)
                ((SamRecordHeaderFilter)filter).setHeader(this.getSAMFileHeader());

        sampleDataSource = new SampleDataSource(getSAMFileHeader(), argCollection.sampleFiles);

        // set the sequence dictionary of all of Tribble tracks to the sequence dictionary of our reference
        rodDataSources = getReferenceOrderedDataSources(referenceMetaDataFiles,referenceDataSource.getReference().getSequenceDictionary(),genomeLocParser,argCollection.unsafe);
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
     * Returns sets of samples present in the (merged) input SAM stream, grouped by readers (i.e. underlying
     * individual bam files). For instance: if GATK is run with three input bam files (three -I arguments), then the list
     * returned by this method will contain 3 elements (one for each reader), with each element being a set of sample names
     * found in the corresponding bam file.
     *
     * @return Sets of samples in the merged input SAM stream, grouped by readers
     */
    public List<Set<String>> getSamplesByReaders() {
        Collection<SAMReaderID> readers = getReadsDataSource().getReaderIDs();

        List<Set<String>> sample_sets = new ArrayList<Set<String>>(readers.size());

        for (SAMReaderID r : readers) {

            Set<String> samples = new HashSet<String>(1);
            sample_sets.add(samples);

            for (SAMReadGroupRecord g : getReadsDataSource().getHeader(r).getReadGroups()) {
                samples.add(g.getSample());
            }
        }

        return sample_sets;

    }

    /**
     * Returns sets of libraries present in the (merged) input SAM stream, grouped by readers (i.e. underlying
     * individual bam files). For instance: if GATK is run with three input bam files (three -I arguments), then the list
     * returned by this method will contain 3 elements (one for each reader), with each element being a set of library names
     * found in the corresponding bam file.
     *
     * @return Sets of libraries present in the (merged) input SAM stream, grouped by readers
     */
    public List<Set<String>> getLibrariesByReaders() {


        Collection<SAMReaderID> readers = getReadsDataSource().getReaderIDs();

        List<Set<String>> lib_sets = new ArrayList<Set<String>>(readers.size());

        for (SAMReaderID r : readers) {

            Set<String> libs = new HashSet<String>(2);
            lib_sets.add(libs);

            for (SAMReadGroupRecord g : getReadsDataSource().getHeader(r).getReadGroups()) {
                libs.add(g.getLibrary());
            }
        }

        return lib_sets;

    }

    /**
     * **** UNLESS YOU HAVE GOOD REASON TO, DO NOT USE THIS METHOD; USE getFileToReadGroupIdMapping() INSTEAD ****
     *
     * Returns sets of (remapped) read groups in input SAM stream, grouped by readers (i.e. underlying
     * individual bam files). For instance: if GATK is run with three input bam files (three -I arguments), then the list
     * returned by this method will contain 3 elements (one for each reader), with each element being a set of remapped read groups
     * (i.e. as seen by read.getReadGroup().getReadGroupId() in the merged stream) that come from the corresponding bam file.
     *
     * @return sets of (merged) read group ids in order of input bams
     */
    public List<Set<String>> getMergedReadGroupsByReaders() {


        Collection<SAMReaderID> readers = getReadsDataSource().getReaderIDs();

        List<Set<String>> rg_sets = new ArrayList<Set<String>>(readers.size());

        for (SAMReaderID r : readers) {

            Set<String> groups = new HashSet<String>(5);
            rg_sets.add(groups);

            for (SAMReadGroupRecord g : getReadsDataSource().getHeader(r).getReadGroups()) {
                if (getReadsDataSource().hasReadGroupCollisions()) { // Check if there were read group clashes with hasGroupIdDuplicates and if so:
                    // use HeaderMerger to translate original read group id from the reader into the read group id in the
                    // merged stream, and save that remapped read group id to associate it with specific reader
                    groups.add(getReadsDataSource().getReadGroupId(r, g.getReadGroupId()));
                } else {
                    // otherwise, pass through the unmapped read groups since this is what Picard does as well
                    groups.add(g.getReadGroupId());
                }
            }
        }

        return rg_sets;

    }

    /**
     * Now that all files are open, validate the sequence dictionaries of the reads vs. the reference vrs the reference ordered data (if available).
     *
     * @param reads     Reads data source.
     * @param reference Reference data source.
     * @param rods    a collection of the reference ordered data tracks
     */
    private void validateSourcesAgainstReference(SAMDataSource reads, ReferenceSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, RMDTrackBuilder manager) {
        if ((reads.isEmpty() && (rods == null || rods.isEmpty())) || reference == null )
            return;

        // Compile a set of sequence names that exist in the reference file.
        SAMSequenceDictionary referenceDictionary = reference.getSequenceDictionary();

        if (!reads.isEmpty()) {
            // Compile a set of sequence names that exist in the BAM files.
            SAMSequenceDictionary readsDictionary = reads.getHeader().getSequenceDictionary();

            Set<String> readsSequenceNames = new TreeSet<String>();
            for (SAMSequenceRecord dictionaryEntry : readsDictionary.getSequences())
                readsSequenceNames.add(dictionaryEntry.getSequenceName());


            if (readsSequenceNames.size() == 0) {
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
     * @return A data source for the given set of reads.
     */
    private SAMDataSource createReadsDataSource(GenomeLocParser genomeLocParser, IndexedFastaSequenceFile refReader) {
        DownsamplingMethod method = getDownsamplingMethod();

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
                argCollection.defaultBaseQualities);
    }

    /**
     * Opens a reference sequence file paired with an index.  Only public for testing purposes
     *
     * @param refFile Handle to a reference sequence file.  Non-null.
     * @return A thread-safe file wrapper.
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
        // try and make the tracks given their requests
        // create of live instances of the tracks
        List<RMDTrack> tracks = new ArrayList<RMDTrack>();

        List<ReferenceOrderedDataSource> dataSources = new ArrayList<ReferenceOrderedDataSource>();
        for (RMDTriplet fileDescriptor : referenceMetaDataFiles)
            dataSources.add(new ReferenceOrderedDataSource(fileDescriptor,
                                                           builder,
                                                           sequenceDictionary,
                                                           genomeLocParser,
                                                           flashbackData()));

        // validation: check to make sure everything the walker needs is present, and that all sequence dictionaries match.
        validateSuppliedReferenceOrderedData(dataSources);
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
     * @return List of intervals.
     */
    public GenomeLocSortedSet getIntervals() {
        return this.intervals;
    }

    /**
     * Gets the list of filters employed by this engine.
     * @return Collection of filters (actual instances) used by this engine.
     */
    public Collection<SamRecordFilter> getFilters() {
        return this.filters;
    }

    /**
     * Sets the list of filters employed by this engine.
     * @param filters Collection of filters (actual instances) used by this engine.
     */
    public void setFilters(Collection<SamRecordFilter> filters) {
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
     * @return cumulative metrics about the entire run.
     */
    public ReadMetrics getCumulativeMetrics() {
        return readsDataSource == null ? null : readsDataSource.getCumulativeReadMetrics();
    }

    public SampleDataSource getSampleMetadata() {
        return this.sampleDataSource;
    }

    /**
     * Get a sample by its ID
     * If an alias is passed in, return the main sample object
     * @param id sample id
     * @return sample Object with this ID
     */
    public Sample getSampleById(String id) {
        return sampleDataSource.getSampleById(id);
    }

    /**
     * Get the sample for a given read group
     * Must first look up ID for read group
     * @param readGroup of sample
     * @return sample object with ID from the read group
     */
    public Sample getSampleByReadGroup(SAMReadGroupRecord readGroup) {
        return sampleDataSource.getSampleByReadGroup(readGroup);
    }

    /**
     * Get a sample for a given read
     * Must first look up read group, and then sample ID for that read group
     * @param read of sample
     * @return sample object of this read
     */
    public Sample getSampleByRead(SAMRecord read) {
        return getSampleByReadGroup(read.getReadGroup());
    }

    /**
     * Get number of sample objects
     * @return size of samples map
     */
    public int sampleCount() {
        return sampleDataSource.sampleCount();
    }

    /**
     * Return all samples with a given family ID
     * Note that this isn't terribly efficient (linear) - it may be worth adding a new family ID data structure for this
     * @param familyId family ID
     * @return Samples with the given family ID
     */
    public Set<Sample> getFamily(String familyId) {
        return sampleDataSource.getFamily(familyId);
    }

    /**
     * Returns all children of a given sample
     * See note on the efficiency of getFamily() - since this depends on getFamily() it's also not efficient
     * @param sample parent sample
     * @return children of the given sample
     */
    public Set<Sample> getChildren(Sample sample) {
        return sampleDataSource.getChildren(sample);
    }

    /**
     * Gets all the samples
     * @return
     */
    public Collection<Sample> getSamples() {
        return sampleDataSource.getSamples();
    }

    /**
     * Takes a list of sample names and returns their corresponding sample objects
     *
     * @param sampleNameList List of sample names
     * @return Corresponding set of samples
     */
    public Set<Sample> getSamples(Collection<String> sampleNameList) {
	return sampleDataSource.getSamples(sampleNameList);
    }


    /**
     * Returns a set of samples that have any value (which could be null) for a given property
     * @param key Property key
     * @return Set of samples with the property
     */
    public Set<Sample> getSamplesWithProperty(String key) {
        return sampleDataSource.getSamplesWithProperty(key);
    }

    /**
     * Returns a set of samples that have a property with a certain value
     * Value must be a string for now - could add a similar method for matching any objects in the future
     *
     * @param key Property key
     * @param value String property value
     * @return Set of samples that match key and value
     */
    public Set<Sample> getSamplesWithProperty(String key, String value) {
        return sampleDataSource.getSamplesWithProperty(key, value);

    }

    /**
     * Returns a set of sample objects for the sample names in a variant context
     *
     * @param context Any variant context
     * @return a set of the sample objects
     */
    public Set<Sample> getSamplesByVariantContext(VariantContext context) {
        Set<Sample> samples = new HashSet<Sample>();
        for (String sampleName : context.getSampleNames()) {
            samples.add(sampleDataSource.getOrCreateSample(sampleName));
        }
        return samples;
    }

    /**
     * Returns all samples that were referenced in the SAM file
     */
    public Set<Sample> getSAMFileSamples() {
        return sampleDataSource.getSAMFileSamples();
    }

    /**
     * Return a subcontext restricted to samples with a given property key/value
     * Gets the sample names from key/value and relies on VariantContext.subContextFromGenotypes for the filtering
     * @param context VariantContext to filter
     * @param key property key
     * @param value property value (must be string)
     * @return subcontext
     */
    public VariantContext subContextFromSampleProperty(VariantContext context, String key, String value) {
        return sampleDataSource.subContextFromSampleProperty(context, key, value);
    }

    public Map<String,String> getApproximateCommandLineArguments(Object... argumentProviders) {
        return CommandLineUtils.getApproximateCommandLineArguments(parsingEngine,argumentProviders);
    }

    public String createApproximateCommandLineArgumentString(Object... argumentProviders) {
        return CommandLineUtils.createApproximateCommandLineArgumentString(parsingEngine,argumentProviders);
    }
    

}
