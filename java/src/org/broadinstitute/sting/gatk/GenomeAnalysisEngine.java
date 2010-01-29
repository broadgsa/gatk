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

package org.broadinstitute.sting.gatk;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.*;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.IndexDrivenSAMDataSource;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.MonolithicShardStrategy;
import org.broadinstitute.sting.gatk.executive.MicroScheduler;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.filters.FilterManager;
import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.bed.BedParser;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.cmdLine.ArgumentException;
import org.broadinstitute.sting.utils.cmdLine.ArgumentSource;
import org.broadinstitute.sting.gatk.io.stubs.Stub;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class GenomeAnalysisEngine {

    // our instance of this genome analysis toolkit; it's used by other classes to extract the traversal engine
    // TODO: public static without final tends to indicate we're thinking about this the wrong way
    public static GenomeAnalysisEngine instance;

    /**
     * Accessor for sharded read data.
     */
    private SAMDataSource readsDataSource = null;

    /**
     * Accessor for sharded reference data.
     */
    private IndexedFastaSequenceFile referenceDataSource = null;

    /**
     * Accessor for sharded reference-ordered data.
     */
    private List<ReferenceOrderedDataSource> rodDataSources;

    // our argument collection
    private GATKArgumentCollection argCollection;

    /**
     * Collection of inputs used by the walker.
     */
    private Map<ArgumentSource, Object> inputs = new HashMap<ArgumentSource, Object>();

    /**
     * Collection of outputs used by the walker.
     */
    private Collection<Stub<?>> outputs = new ArrayList<Stub<?>>();

    /**
     * Collection of the filters applied to the walker's input data.
     */
    private Collection<SamRecordFilter> filters;

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(GenomeAnalysisEngine.class);

    /**
     * our walker manager
     */
    private final WalkerManager walkerManager;

    /**
     * Manage lists of filters.
     */
    private final FilterManager filterManager;

    /**
     * our constructor, where all the work is done
     * <p/>
     * legacy traversal types are sent to legacyTraversal function; as we move more of the traversals to the
     * new MicroScheduler class we'll be able to delete that function.
     */
    public GenomeAnalysisEngine() {
        // make sure our instance variable points to this analysis engine
        instance = this;
        walkerManager = new WalkerManager();
        filterManager = new FilterManager();
    }

    /**
     * Actually run the GATK with the specified walker.
     *
     * @param args      the argument collection, where we get all our setup information from
     * @param my_walker Walker to run over the dataset.  Must not be null.
     * @return the value of this traversal.
     */
    public Object execute(GATKArgumentCollection args, Walker<?, ?> my_walker, Collection<SamRecordFilter> filters) {
        // validate our parameters
        if (args == null) {
            throw new StingException("The GATKArgumentCollection passed to GenomeAnalysisEngine can not be null.");
        }

        // validate our parameters
        if (my_walker == null)
            throw new StingException("The walker passed to GenomeAnalysisEngine can not be null.");

        // save our argument parameter
        this.argCollection = args;
        this.filters = filters;

        // Prepare the data for traversal.
        initializeDataSources(my_walker, filters, argCollection);

        // our microscheduler, which is in charge of running everything
        MicroScheduler microScheduler = createMicroscheduler(my_walker);

        // create the output streams
        initializeOutputStreams(my_walker, microScheduler.getOutputTracker());

        GenomeLocSortedSet locs = null;
        if (argCollection.intervals != null && argCollection.intervalMerging.check()) {
            locs = GenomeLocSortedSet.createSetFromList(parseIntervalRegion(argCollection.intervals));
        }

        ShardStrategy shardStrategy = getShardStrategy(my_walker,
                                                       microScheduler.getReference(),
                                                       locs,
                                                       argCollection.maximumEngineIterations,
                                                       readsDataSource != null ? readsDataSource.getReadsInfo().getValidationExclusionList() : null);

        // execute the microscheduler, storing the results
        return microScheduler.execute(my_walker, shardStrategy, argCollection.maximumEngineIterations);
    }

    /**
     * Add additional, externally managed IO streams for walker input.
     *
     * @param argumentSource Field in the walker into which to inject the value.
     * @param value          Instance to inject.
     */
    public void addInput(ArgumentSource argumentSource, Object value) {
        inputs.put(argumentSource, value);
    }

    /**
     * Add additional, externally managed IO streams for walker output.
     *
     * @param stub Instance to inject.
     */
    public void addOutput(Stub<?> stub) {
        outputs.add(stub);
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
    public String getWalkerName(Class<Walker> walkerType) {
        return walkerManager.getName(walkerType);
    }

    /**
     * Gets a list of the filters to associate with the given walker.  Will NOT initialize the engine with this filters;
     * the caller must handle that directly.
     * @param args Existing argument collection, for compatibility with legacy command-line walkers.
     * @param walker Walker to use when determining which filters to apply.
     * @return A collection of available filters.
     */
    protected Collection<SamRecordFilter> createFiltersForWalker(GATKArgumentCollection args, Walker walker) {
        Set<SamRecordFilter> filters = new HashSet<SamRecordFilter>();        
        filters.addAll(WalkerManager.getReadFilters(walker,filterManager));
        if (args.filterZeroMappingQualityReads != null && args.filterZeroMappingQualityReads)
            filters.add(new ZeroMappingQualityReadFilter());
        for(String filterName: args.readFilters)
            filters.add(filterManager.createByName(filterName));
        return Collections.unmodifiableSet(filters);
    }

    /**
     * Allow subclasses and others within this package direct access to the walker manager.
     * @return The walker manager used by this package.
     */
    protected WalkerManager getWalkerManager() {
        return walkerManager;
    }

    private void initializeDataSources(Walker my_walker, Collection<SamRecordFilter> filters, GATKArgumentCollection argCollection) {
        validateSuppliedReadsAgainstWalker(my_walker, argCollection);
        logger.info("Strictness is " + argCollection.strictnessLevel);
        readsDataSource = createReadsDataSource(extractSourceInfo(my_walker, filters, argCollection));

        validateSuppliedReferenceAgainstWalker(my_walker, argCollection);
        referenceDataSource = openReferenceSequenceFile(argCollection.referenceFile);

        validateReadsAndReferenceAreCompatible(readsDataSource, referenceDataSource);

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
        // TODO: The ROD iterator currently does not understand multiple intervals file.  Fix this by cleaning the ROD system.
        if (argCollection.intervals != null && argCollection.intervals.size() == 1) {
            bindConvenienceRods("interval", "Intervals", argCollection.intervals.get(0).replaceAll(",", ""));
        }

        // parse out the rod bindings
        ReferenceOrderedData.parseBindings(argCollection.RODBindings, rods);

        validateSuppliedReferenceOrderedDataAgainstWalker(my_walker, rods);

        rodDataSources = getReferenceOrderedDataSources(rods);
    }

    /**
     * setup a microscheduler
     *
     * @param my_walker our walker of type LocusWalker
     * @return a new microscheduler
     */
    private MicroScheduler createMicroscheduler(Walker my_walker) {
        // the mircoscheduler to return
        MicroScheduler microScheduler = null;

        // we need to verify different parameter based on the walker type
        if (my_walker instanceof LocusWalker || my_walker instanceof LocusWindowWalker) {
            // create the MicroScheduler
            microScheduler = MicroScheduler.create(my_walker, readsDataSource, referenceDataSource, rodDataSources, argCollection.numberOfThreads);
        } else if (my_walker instanceof ReadWalker || my_walker instanceof DuplicateWalker) {
            if (argCollection.referenceFile == null)
                Utils.scareUser(String.format("Read-based traversals require a reference file but none was given"));
            microScheduler = MicroScheduler.create(my_walker, readsDataSource, referenceDataSource, rodDataSources, argCollection.numberOfThreads);
        } else {
            Utils.scareUser(String.format("Unable to create the appropriate TraversalEngine for analysis type %s", walkerManager.getName(my_walker.getClass())));
        }

        return microScheduler;
    }

    /**
     * setup the interval regions, from either the interval file of the genome region string
     *
     * @param intervals the list of intervals to parse
     * @return a list of genomeLoc representing the interval file
     */
    public static List<GenomeLoc> parseIntervalRegion(final List<String> intervals) {
        List<GenomeLoc> locs = new ArrayList<GenomeLoc>();
        for (String interval : intervals) {
            if (new File(interval).exists()) {
                // support for the bed style interval format
                if (interval.toUpperCase().endsWith(".BED")) {
                    Utils.warnUser("Bed files are 0 based half-open intervals, which are converted to 1-based closed intervals in the GATK.  " +
                            "Be aware that all output information and intervals are 1-based closed intervals.");
                    BedParser parser = new BedParser(new File(interval));
                    locs.addAll(parser.getSortedAndMergedLocations(GenomeAnalysisEngine.instance.getArguments().intervalMerging));
                } else {
                    locs.addAll(GenomeLocParser.intervalFileToList(interval,GenomeAnalysisEngine.instance.getArguments().intervalMerging));
                }
            } else {
                locs.addAll(GenomeLocParser.parseGenomeLocs(interval,GenomeAnalysisEngine.instance.getArguments().intervalMerging));
            }

        }
        return locs;
    }

    /**
     * Returns sets of samples present in the (merged) input SAM stream, grouped by readers (i.e. underlying
     * individual bam files). For instance: if GATK is run with three input bam files (three -I arguments), then the list
     * returned by this method will contain 3 elements (one for each reader), with each element being a set of sample names
     * found in the corresponding bam file.
     *
     * @return
     */
    public List<Set<String>> getSamplesByReaders() {


        Collection<SAMFileReader> readers = getDataSource().getReaders();

        List<Set<String>> sample_sets = new ArrayList<Set<String>>(readers.size());

        for (SAMFileReader r : readers) {

            Set<String> samples = new HashSet<String>(1);
            sample_sets.add(samples);

            for (SAMReadGroupRecord g : r.getFileHeader().getReadGroups()) {
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
     * @return
     */
    public List<Set<String>> getLibrariesByReaders() {


        Collection<SAMFileReader> readers = getDataSource().getReaders();

        List<Set<String>> lib_sets = new ArrayList<Set<String>>(readers.size());

        for (SAMFileReader r : readers) {

            Set<String> libs = new HashSet<String>(2);
            lib_sets.add(libs);

            for (SAMReadGroupRecord g : r.getFileHeader().getReadGroups()) {
                libs.add(g.getLibrary());
            }
        }

        return lib_sets;

    }

    /**
     * Returns sets of (remapped) read groups in input SAM stream, grouped by readers (i.e. underlying
     * individual bam files). For instance: if GATK is run with three input bam files (three -I arguments), then the list
     * returned by this method will contain 3 elements (one for each reader), with each element being a set of remapped read groups
     * (i.e. as seen by read.getReadGroup().getReadGroupId() in the merged stream) that come from the corresponding bam file.
     *
     * @return
     */
    public List<Set<String>> getMergedReadGroupsByReaders() {


        Collection<SAMFileReader> readers = getDataSource().getReaders();

        List<Set<String>> rg_sets = new ArrayList<Set<String>>(readers.size());

        for (SAMFileReader r : readers) {

            Set<String> groups = new HashSet<String>(5);
            rg_sets.add(groups);

            for (SAMReadGroupRecord g : r.getFileHeader().getReadGroups()) {
                if (getDataSource().hasReadGroupCollisions()) { // Check if there were read group clashes with hasGroupIdDuplicates and if so:
                    // use HeaderMerger to translate original read group id from the reader into the read group id in the
                    // merged stream, and save that remapped read group id to associate it with specific reader
                    groups.add(getDataSource().getReadGroupId(r, g.getReadGroupId()));
                } else {
                    // otherwise, pass through the unmapped read groups since this is what Picard does as well
                    groups.add(g.getReadGroupId());
                }
            }
        }

        return rg_sets;

    }

    /**
     * Bundles all the source information about the reads into a unified data structure.
     *
     * @param walker        The walker for which to extract info.
     * @param argCollection The collection of arguments passed to the engine.
     * @return The reads object providing reads source info.
     */
    private Reads extractSourceInfo(Walker walker, Collection<SamRecordFilter> filters, GATKArgumentCollection argCollection) {
        return new Reads(argCollection.samFiles,
                argCollection.strictnessLevel,
                argCollection.downsampleFraction,
                argCollection.downsampleCoverage,
                new ValidationExclusion(Arrays.asList(argCollection.unsafe)),
                filters,
                argCollection.readMaxPileup,
                walker.includeReadsWithDeletionAtLoci(),
                walker.generateExtendedEvents());
    }

    /**
     * Verifies that the supplied set of reads files mesh with what the walker says it requires.
     *
     * @param walker    Walker to test.
     * @param arguments Supplied reads files.
     */
    private void validateSuppliedReadsAgainstWalker(Walker walker, GATKArgumentCollection arguments) {
        // Check what the walker says is required against what was provided on the command line.
        if (WalkerManager.isRequired(walker, DataSource.READS) && (arguments.samFiles == null || arguments.samFiles.size() == 0))
            throw new ArgumentException("Walker requires reads but none were provided.  If this is incorrect, alter the walker's @Requires annotation.");

        // Check what the walker says is allowed against what was provided on the command line.
        if ((arguments.samFiles != null && arguments.samFiles.size() > 0) && !WalkerManager.isAllowed(walker, DataSource.READS))
            throw new ArgumentException("Walker does not allow reads but reads were provided.  If this is incorrect, alter the walker's @Allows annotation");
    }

    /**
     * Verifies that the supplied reference file mesh with what the walker says it requires.
     *
     * @param walker    Walker to test.
     * @param arguments Supplied reads files.
     */
    private void validateSuppliedReferenceAgainstWalker(Walker walker, GATKArgumentCollection arguments) {
        // Check what the walker says is required against what was provided on the command line.
        // TODO: Temporarily disabling WalkerManager.isRequired check on the reference because the reference is always required.
        if (/*WalkerManager.isRequired(walker, DataSource.REFERENCE) &&*/ arguments.referenceFile == null)
            throw new ArgumentException("Walker requires a reference but none was provided.  If this is incorrect, alter the walker's @Requires annotation.");

        // Check what the walker says is allowed against what was provided on the command line.
        if (arguments.referenceFile != null && !WalkerManager.isAllowed(walker, DataSource.REFERENCE))
            throw new ArgumentException("Walker does not allow a reference but one was provided.  If this is incorrect, alter the walker's @Allows annotation");
    }

    /**
     * Verifies that all required reference-ordered data has been supplied, and any reference-ordered data that was not
     * 'allowed' is still present.
     *
     * @param walker Walker to test.
     * @param rods   Reference-ordered data to load.
     */
    private void validateSuppliedReferenceOrderedDataAgainstWalker(Walker walker, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        // Check to make sure that all required metadata is present.
        List<RMD> allRequired = WalkerManager.getRequiredMetaData(walker);
        for (RMD required : allRequired) {
            boolean found = false;
            for (ReferenceOrderedData<? extends ReferenceOrderedDatum> rod : rods) {
                if (rod.matches(required.name(), required.type()))
                    found = true;
            }
            if (!found)
                throw new ArgumentException(String.format("Unable to find reference metadata (%s,%s)", required.name(), required.type()));
        }

        // Check to see that no forbidden rods are present.
        for (ReferenceOrderedData<? extends ReferenceOrderedDatum> rod : rods) {
            if (!WalkerManager.isAllowed(walker, rod))
                throw new ArgumentException(String.format("Walker of type %s does not allow access to metadata: %s.  If this is incorrect, change the @Allows metadata", walker.getClass(), rod.getName()));
        }
    }

    /**
     * Now that all files are open, validate the sequence dictionaries of the reads vs. the reference.
     *
     * @param reads     Reads data source.
     * @param reference Reference data source.
     */
    private void validateReadsAndReferenceAreCompatible(SAMDataSource reads, ReferenceSequenceFile reference) {
        if (reads == null || reference == null)
            return;

        // Compile a set of sequence names that exist in the BAM files.
        SAMSequenceDictionary readsDictionary = reads.getHeader().getSequenceDictionary();

        Set<String> readsSequenceNames = new TreeSet<String>();
        for (SAMSequenceRecord dictionaryEntry : readsDictionary.getSequences())
            readsSequenceNames.add(dictionaryEntry.getSequenceName());

        // Compile a set of sequence names that exist in the reference file.
        SAMSequenceDictionary referenceDictionary = reference.getSequenceDictionary();

        Set<String> referenceSequenceNames = new TreeSet<String>();
        for (SAMSequenceRecord dictionaryEntry : referenceDictionary.getSequences())
            referenceSequenceNames.add(dictionaryEntry.getSequenceName());

        if (readsSequenceNames.size() == 0) {
            logger.info("Reads file is unmapped.  Skipping validation against reference.");
            return;
        }

        // If there's no overlap between reads and reference, data will be bogus.  Throw an exception.
        Set<String> intersectingSequenceNames = new HashSet<String>(readsSequenceNames);
        intersectingSequenceNames.retainAll(referenceSequenceNames);
        if (intersectingSequenceNames.size() == 0) {
            StringBuilder error = new StringBuilder();
            error.append("No overlap exists between sequence dictionary of the reads and the sequence dictionary of the reference.  Perhaps you're using the wrong reference?\n");
            error.append(System.getProperty("line.separator"));
            error.append(String.format("Reads contigs:     %s%n", prettyPrintSequenceRecords(readsDictionary)));
            error.append(String.format("Reference contigs: %s%n", prettyPrintSequenceRecords(referenceDictionary)));
            logger.error(error.toString());
            Utils.scareUser("No overlap exists between sequence dictionary of the reads and the sequence dictionary of the reference.");
        }

        // If the two datasets are not equal and neither is a strict subset of the other, warn the user.
        if (!readsSequenceNames.equals(referenceSequenceNames) &&
                !readsSequenceNames.containsAll(referenceSequenceNames) &&
                !referenceSequenceNames.containsAll(readsSequenceNames)) {
            StringBuilder warning = new StringBuilder();
            warning.append("Limited overlap exists between sequence dictionary of the reads and the sequence dictionary of the reference.  Perhaps you're using the wrong reference?\n");
            warning.append(System.getProperty("line.separator"));
            warning.append(String.format("Reads contigs:     %s%n", prettyPrintSequenceRecords(readsDictionary)));
            warning.append(String.format("Reference contigs: %s%n", prettyPrintSequenceRecords(referenceDictionary)));
            logger.warn(warning.toString());
        }
    }

    private String prettyPrintSequenceRecords(SAMSequenceDictionary sequenceDictionary) {
        String[] sequenceRecordNames = new String[sequenceDictionary.size()];
        int sequenceRecordIndex = 0;
        for (SAMSequenceRecord sequenceRecord : sequenceDictionary.getSequences())
            sequenceRecordNames[sequenceRecordIndex++] = sequenceRecord.getSequenceName();
        return Arrays.deepToString(sequenceRecordNames);
    }


    /**
     * Convenience function that binds RODs using the old-style command line parser to the new style list for
     * a uniform processing.
     *
     * @param name the name of the rod
     * @param type its type
     * @param file the file to load the rod from
     */
    private void bindConvenienceRods(final String name, final String type, final String file) {
        argCollection.RODBindings.add(Utils.join(",", new String[]{name, type, file}));
    }

    /**
     * Get the sharding strategy given a driving data source.
     *
     * @param walker            Walker for which to infer sharding strategy.
     * @param drivingDataSource Data on which to shard.
     * @param intervals         Intervals to use when limiting sharding.
     * @param maxIterations     the maximum number of iterations to run through
     * @return Sharding strategy for this driving data source.
     */
    protected ShardStrategy getShardStrategy(Walker walker,
                                             ReferenceSequenceFile drivingDataSource,
                                             GenomeLocSortedSet intervals,
                                             Integer maxIterations,
                                             ValidationExclusion exclusions) {
        if(readsDataSource != null && !readsDataSource.hasIndex()) {
            if(!exclusions.contains(ValidationExclusion.TYPE.ALLOW_UNINDEXED_BAM) || intervals != null)
                throw new StingException("The GATK cannot currently process unindexed BAM files");

            Shard.ShardType shardType;
            if(walker instanceof LocusWalker)
                shardType = Shard.ShardType.LOCUS;
            else if(walker instanceof ReadWalker || walker instanceof DuplicateWalker)
                shardType = Shard.ShardType.READ;
            else
                throw new StingException("The GATK cannot currently process unindexed BAM files");

            return new MonolithicShardStrategy(shardType);
        }

        ShardStrategy shardStrategy = null;
        ShardStrategyFactory.SHATTER_STRATEGY shardType;

        long SHARD_SIZE = 100000L;

        if (walker instanceof LocusWalker) {
            if (walker instanceof RodWalker) SHARD_SIZE *= 1000;

            if (intervals != null && !intervals.isEmpty()) {
                shardType = (walker.isReduceByInterval()) ?
                        ShardStrategyFactory.SHATTER_STRATEGY.INTERVAL :
                        ShardStrategyFactory.SHATTER_STRATEGY.LINEAR;
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                        argCollection.experimentalSharding ? ShardStrategyFactory.SHATTER_STRATEGY.LOCUS_EXPERIMENTAL : shardType,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE,
                        intervals, maxIterations);
            } else
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,ShardStrategyFactory.SHATTER_STRATEGY.LINEAR,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE, maxIterations);
        } else if (walker instanceof ReadWalker ||
                walker instanceof DuplicateWalker) {
            if(argCollection.experimentalSharding)
                shardType = ShardStrategyFactory.SHATTER_STRATEGY.READS_EXPERIMENTAL;
            else
                shardType = ShardStrategyFactory.SHATTER_STRATEGY.READS;

            if (intervals != null && !intervals.isEmpty()) {
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,shardType,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE,
                        intervals, maxIterations);
            } else {
                shardStrategy = ShardStrategyFactory.shatter(readsDataSource,shardType,
                        drivingDataSource.getSequenceDictionary(),
                        SHARD_SIZE, maxIterations);
            }
        } else if (walker instanceof LocusWindowWalker) {
            if ((intervals == null || intervals.isEmpty()) && !exclusions.contains(ValidationExclusion.TYPE.ALLOW_EMPTY_INTERVAL_LIST))
                Utils.warnUser("walker is of type LocusWindow (which operates over intervals), but no intervals were provided." +
                               "This may be unintentional, check your command-line arguments.");
            shardStrategy = ShardStrategyFactory.shatter(readsDataSource,
                    ShardStrategyFactory.SHATTER_STRATEGY.INTERVAL,
                    drivingDataSource.getSequenceDictionary(),
                    SHARD_SIZE,
                    intervals, maxIterations);
        } else
            throw new StingException("Unable to support walker of type" + walker.getClass().getName());

        return shardStrategy;
    }

    /**
     * Gets a data source for the given set of reads.
     *
     * @param reads the read source information
     * @return A data source for the given set of reads.
     */
    private SAMDataSource createReadsDataSource(Reads reads) {
        // By reference traversals are happy with no reads.  Make sure that case is handled.
        if (reads.getReadsFiles().size() == 0)
            return null;

        SAMDataSource dataSource = null;
        if(argCollection.experimentalSharding)
            dataSource = new BlockDrivenSAMDataSource(reads);
        else
            dataSource = new IndexDrivenSAMDataSource(reads);

        return dataSource;
    }

    /**
     * Opens a reference sequence file paired with an index.
     *
     * @param refFile Handle to a reference sequence file.  Non-null.
     * @return A thread-safe file wrapper.
     */
    private IndexedFastaSequenceFile openReferenceSequenceFile(File refFile) {
        IndexedFastaSequenceFile ref = null;
        try {
            ref = new IndexedFastaSequenceFile(refFile);
        }
        catch (FileNotFoundException ex) {
            throw new StingException("I/O error while opening fasta file: " + ex.getMessage(), ex);
        }
        GenomeLocParser.setupRefContigOrdering(ref);
        return ref;
    }

    /**
     * Open the reference-ordered data sources.
     *
     * @param rods the reference order data to execute using
     * @return A list of reference-ordered data sources.
     */
    private List<ReferenceOrderedDataSource> getReferenceOrderedDataSources(List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        List<ReferenceOrderedDataSource> dataSources = new ArrayList<ReferenceOrderedDataSource>();
        for (ReferenceOrderedData<? extends ReferenceOrderedDatum> rod : rods)
            dataSources.add(new ReferenceOrderedDataSource(rod));
        return dataSources;
    }

    /**
     * Initialize the output streams as specified by the user.
     *
     * @param walker        the walker to initialize output streams for
     * @param outputTracker the tracker supplying the initialization data.
     */
    private void initializeOutputStreams(Walker walker, OutputTracker outputTracker) {
        if (argCollection.outErrFileName != null)
            outputTracker.initializeCoreIO(argCollection.outErrFileName, argCollection.outErrFileName);
        else
            outputTracker.initializeCoreIO(argCollection.outFileName, argCollection.errFileName);

        for (Map.Entry<ArgumentSource, Object> input : inputs.entrySet())
            outputTracker.addInput(input.getKey(), input.getValue());
        for (Stub<?> stub : outputs)
            outputTracker.addOutput(stub);

        outputTracker.prepareWalker(walker);
    }

    /**
     * Returns the SAM File Header from the input reads' data source file
     * @return the SAM File Header from the input reads' data source file
     */
    public SAMFileHeader getSAMFileHeader() {
        return readsDataSource.getHeader();
    }

    /**
     * Returns data source object encapsulating all essential info and handlers used to traverse
     * reads; header merger, individual file readers etc can be accessed through the returned data source object.
     *
     * @return the reads data source
     */
    public SAMDataSource getDataSource() {
        return this.readsDataSource;
    }

    /**
     * Gets the collection of GATK main application arguments for enhanced walker validation.
     *
     * @return the GATK argument collection
     */
    public GATKArgumentCollection getArguments() {
        return this.argCollection;
    }

    /**
     * Gets the list of filters employed by this walker.
     * @return Collection of filters (actual instances) used by this walker.
     */
    public Collection<SamRecordFilter> getFilters() {
        return this.filters;
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
}
