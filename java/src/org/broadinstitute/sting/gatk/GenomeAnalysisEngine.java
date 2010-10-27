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
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.commandline.ArgumentException;
import org.broadinstitute.sting.commandline.ArgumentSource;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.shards.MonolithicShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.executive.MicroScheduler;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.io.stubs.Stub;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.*;

/**
 * A GenomeAnalysisEngine that runs a specified walker.
 */
public class GenomeAnalysisEngine extends AbstractGenomeAnalysisEngine {
    /**
     * our walker manager
     */
    private final WalkerManager walkerManager = new WalkerManager();

    private Walker<?, ?> walker;

    public void setWalker(Walker<?, ?> walker) {
        this.walker = walker;
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

        // our microscheduler, which is in charge of running everything
        MicroScheduler microScheduler = createMicroscheduler();

        // create the output streams                     "
        initializeOutputStreams(microScheduler.getOutputTracker());

        initializeIntervals();

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

    /**
     * Gets a list of the filters to associate with the given walker.  Will NOT initialize the engine with this filters;
     * the caller must handle that directly.
     * @return A collection of available filters.
     */
    @Override
    public Collection<SamRecordFilter> createFilters() {
        Set<SamRecordFilter> filters = new HashSet<SamRecordFilter>();
        filters.addAll(WalkerManager.getReadFilters(walker,this.getFilterManager()));
        filters.addAll(super.createFilters());
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
        Walker my_walker = this.walker;

        // the mircoscheduler to return
        MicroScheduler microScheduler = null;

        // Temporarily require all walkers to have a reference, even if that reference is not conceptually necessary.
        if ((my_walker instanceof ReadWalker || my_walker instanceof DuplicateWalker || my_walker instanceof ReadPairWalker) && 
                this.getArguments().referenceFile == null) {
            throw new UserException.CommandLineException("Read-based traversals require a reference file but none was given");
        }

        return MicroScheduler.create(this,my_walker,this.getDataSource(),this.getReferenceDataSource().getReference(),this.getRodDataSources(),this.getArguments().numberOfThreads);
    }

    @Override
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

    @Override
    protected boolean generateExtendedEvents() {
        return walker.generateExtendedEvents();
    }

    @Override
    protected boolean includeReadsWithDeletionAtLoci() {
        return walker.includeReadsWithDeletionAtLoci();
    }

    /**
     * Verifies that the supplied set of reads files mesh with what the walker says it requires.
     */
    @Override
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
    @Override
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
    @Override
    protected void validateSuppliedReferenceOrderedData(List<RMDTrack> rods) {
        // Check to make sure that all required metadata is present.
        List<RMD> allRequired = WalkerManager.getRequiredMetaData(walker);
        for (RMD required : allRequired) {
            boolean found = false;
            for (RMDTrack rod : rods) {
                if (rod.matchesNameAndRecordType(required.name(), required.type()))
                    found = true;
            }
            if (!found)
                throw new ArgumentException(String.format("Walker requires reference metadata to be supplied named '%s' of type '%s', but this metadata was not provided.  " +
                                                          "Please supply the specified metadata file.", required.name(), required.type().getSimpleName()));
        }

        // Check to see that no forbidden rods are present.
        for (RMDTrack rod : rods) {
            if (!WalkerManager.isAllowed(walker, rod))
                throw new ArgumentException(String.format("Walker of type %s does not allow access to metadata: %s", walker.getClass(), rod.getName()));
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
        SAMDataSource readsDataSource = this.getDataSource();
        ValidationExclusion exclusions = (readsDataSource != null ? readsDataSource.getReadsInfo().getValidationExclusionList() : null);
        ReferenceDataSource referenceDataSource = this.getReferenceDataSource();
        // Use monolithic sharding if no index is present.  Monolithic sharding is always required for the original
        // sharding system; it's required with the new sharding system only for locus walkers.
        if(readsDataSource != null && !readsDataSource.hasIndex() ) { 
            if(!exclusions.contains(ValidationExclusion.TYPE.ALLOW_UNINDEXED_BAM))
                throw new UserException.CommandLineException("The GATK cannot currently process unindexed BAM files without the -U ALLOW_UNINDEXED_BAM");
            if(intervals != null && WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE)
                throw new UserException.CommandLineException("Cannot perform interval processing when walker is not driven by reference and no index is available.");

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
                    region.add(GenomeLocParser.createGenomeLoc(sequenceRecord.getSequenceName(),1,sequenceRecord.getSequenceLength()));
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
                        intervals);
            } else
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        referenceDataSource.getReference(),
                        ShardStrategyFactory.SHATTER_STRATEGY.LOCUS_EXPERIMENTAL,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE);
        } else if (walker instanceof ReadWalker ||
                walker instanceof DuplicateWalker) {
            shardType = ShardStrategyFactory.SHATTER_STRATEGY.READS_EXPERIMENTAL;

            if (intervals != null && !intervals.isEmpty()) {
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        referenceDataSource.getReference(),
                        shardType,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE,
                        intervals);
            } else {
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        referenceDataSource.getReference(),
                        shardType,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE);
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
                    SHARD_SIZE);
        } else
            throw new ReviewedStingException("Unable to support walker of type" + walker.getClass().getName());

        return shardStrategy;
    }

    @Override
    protected boolean flashbackData() {
        return walker instanceof ReadWalker;
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

}
