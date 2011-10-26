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

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.DownsamplingMethod;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.filters.CountingFilteringIterator;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.baq.BAQSamIterator;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSamRecordFactory;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 2:36:16 PM
 * <p/>
 * Converts shards to SAM iterators over the specified region
 */
public class SAMDataSource {
    final private static GATKSamRecordFactory factory = new GATKSamRecordFactory();

    /** Backing support for reads. */
    protected final ReadProperties readProperties;

    /**
     * Runtime metrics of reads filtered, etc.
     */
    private final ReadMetrics readMetrics;

    /**
     * Tools for parsing GenomeLocs, for verifying BAM ordering against general ordering.
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * Identifiers for the readers driving this data source.
     */
    private final Collection<SAMReaderID> readerIDs;

    /**
     * How strict are the readers driving this data source.
     */
    private final SAMFileReader.ValidationStringency validationStringency;

    /**
     * Store BAM indices for each reader present.
     */
    private final Map<SAMReaderID,GATKBAMIndex> bamIndices = new HashMap<SAMReaderID,GATKBAMIndex>();

    /**
     * How far along is each reader?
     */
    private final Map<SAMReaderID, SAMFileSpan> readerPositions = new HashMap<SAMReaderID,SAMFileSpan>();

    /**
     * The merged header.
     */
    private final SAMFileHeader mergedHeader;

    /**
     * The sort order of the BAM files.  Files without a sort order tag are assumed to be
     * in coordinate order.
     */
    private SAMFileHeader.SortOrder sortOrder = null;

    /**
     * Whether the read groups in overlapping files collide.
     */
    private final boolean hasReadGroupCollisions;

    /**
     * Maps the SAM readers' merged read group ids to their original ids. Since merged read group ids
     * are always unique, we can simply use a map here, no need to stratify by reader.
     */
    private final ReadGroupMapping mergedToOriginalReadGroupMappings = new ReadGroupMapping();

    /**
     * Maps the SAM readers' original read group ids to their revised ids. This mapping must be stratified
     * by readers, since there can be readgroup id collision: different bam files (readers) can list the
     * same read group id, which will be disambiguated when these input streams are merged.
     */
    private final Map<SAMReaderID,ReadGroupMapping> originalToMergedReadGroupMappings = new HashMap<SAMReaderID,ReadGroupMapping>();

    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(SAMDataSource.class);

    /**
     * A collection of readers driving the merging process.
     */
    private final SAMResourcePool resourcePool;

    /**
     * Whether to enable the new low-memory sharding mechanism.
     */
    private boolean enableLowMemorySharding = false;

    /**
     * Create a new SAM data source given the supplied read metadata.
     * @param samFiles list of reads files.
     */
    public SAMDataSource(Collection<SAMReaderID> samFiles,GenomeLocParser genomeLocParser) {
        this(
                samFiles,
                genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                false,
                false,
                true);
    }

    /**
     * See complete constructor.  Does not enable BAQ by default.
     */
    public SAMDataSource(
            Collection<SAMReaderID> samFiles,
            GenomeLocParser genomeLocParser,
            boolean useOriginalBaseQualities,
            SAMFileReader.ValidationStringency strictness,
            Integer readBufferSize,
            DownsamplingMethod downsamplingMethod,
            ValidationExclusion exclusionList,
            Collection<ReadFilter> supplementalFilters,
            boolean includeReadsWithDeletionAtLoci,
            boolean generateExtendedEvents,
            boolean enableLowMemorySharding) {
        this(   samFiles,
                genomeLocParser,
                useOriginalBaseQualities,
                strictness,
                readBufferSize,
                downsamplingMethod,
                exclusionList,
                supplementalFilters,
                includeReadsWithDeletionAtLoci,
                generateExtendedEvents,
                BAQ.CalculationMode.OFF,
                BAQ.QualityMode.DONT_MODIFY,
                null, // no BAQ
                (byte) -1,
                enableLowMemorySharding);
        }

    /**
     * Create a new SAM data source given the supplied read metadata.
     * @param samFiles list of reads files.
     * @param useOriginalBaseQualities True if original base qualities should be used.
     * @param strictness Stringency of reads file parsing.
     * @param readBufferSize Number of reads to hold in memory per BAM.
     * @param downsamplingMethod Method for downsampling reads at a given locus.
     * @param exclusionList what safety checks we're willing to let slide
     * @param supplementalFilters additional filters to dynamically apply.
     * @param generateExtendedEvents if true, the engine will issue an extra call to walker's map() with
     *        a pile of indel/noevent extended events at every locus with at least one indel associated with it
     *        (in addition to a "regular" call to map() at this locus performed with base pileup)
     * @param includeReadsWithDeletionAtLoci if 'true', the base pileups sent to the walker's map() method
     *         will explicitly list reads with deletion over the current reference base; otherwise, only observed
     *        bases will be seen in the pileups, and the deletions will be skipped silently.
     * @param defaultBaseQualities if the reads have incomplete quality scores, set them all to defaultBaseQuality.
     */
    public SAMDataSource(
            Collection<SAMReaderID> samFiles,
            GenomeLocParser genomeLocParser,
            boolean useOriginalBaseQualities,
            SAMFileReader.ValidationStringency strictness,
            Integer readBufferSize,
            DownsamplingMethod downsamplingMethod,
            ValidationExclusion exclusionList,
            Collection<ReadFilter> supplementalFilters,
            boolean includeReadsWithDeletionAtLoci,
            boolean generateExtendedEvents,
            BAQ.CalculationMode cmode,
            BAQ.QualityMode qmode,
            IndexedFastaSequenceFile refReader,
            byte defaultBaseQualities,
            boolean enableLowMemorySharding) {
        this.enableLowMemorySharding(enableLowMemorySharding);
        this.readMetrics = new ReadMetrics();
        this.genomeLocParser = genomeLocParser;

        readerIDs = samFiles;
        validationStringency = strictness;
        for (SAMReaderID readerID : samFiles) {
            if (!readerID.samFile.canRead())
                throw new UserException.CouldNotReadInputFile(readerID.samFile,"file is not present or user does not have appropriate permissions.  " +
                                                                               "Please check that the file is present and readable and try again.");
        }

        resourcePool = new SAMResourcePool(Integer.MAX_VALUE);
        SAMReaders readers = resourcePool.getAvailableReaders();

        // Determine the sort order.
        for(SAMFileReader reader: readers.values()) {
            // Get the sort order, forcing it to coordinate if unsorted.
            SAMFileHeader header = reader.getFileHeader();

            if ( header.getReadGroups().isEmpty() ) {
                throw new UserException.MalformedBAM(readers.getReaderID(reader).samFile,
                        "SAM file doesn't have any read groups defined in the header.  The GATK no longer supports SAM files without read groups");
            }

            SAMFileHeader.SortOrder sortOrder = header.getSortOrder() != SAMFileHeader.SortOrder.unsorted ? header.getSortOrder() : SAMFileHeader.SortOrder.coordinate;

            // Validate that all input files are sorted in the same order.
            if(this.sortOrder != null && this.sortOrder != sortOrder)
                throw new UserException.MissortedBAM(String.format("Attempted to process mixed of files sorted as %s and %s.",this.sortOrder,sortOrder));

            // Update the sort order.
            this.sortOrder = sortOrder;
        }

        initializeReaderPositions(readers);

        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate,readers.headers(),true);
        mergedHeader = headerMerger.getMergedHeader();
        hasReadGroupCollisions = headerMerger.hasReadGroupCollisions();

        readProperties = new ReadProperties(
                samFiles,
                mergedHeader,
                useOriginalBaseQualities,
                strictness,
                readBufferSize,
                downsamplingMethod,
                exclusionList,
                supplementalFilters,
                includeReadsWithDeletionAtLoci,
                generateExtendedEvents,
                cmode,
                qmode,
                refReader,
                defaultBaseQualities);
        
        // cache the read group id (original) -> read group id (merged)
        // and read group id (merged) -> read group id (original) mappings.
        for(SAMReaderID id: readerIDs) {
            SAMFileReader reader = readers.getReader(id);
            ReadGroupMapping mappingToMerged = new ReadGroupMapping();

            List<SAMReadGroupRecord> readGroups = reader.getFileHeader().getReadGroups();
            for(SAMReadGroupRecord readGroup: readGroups) {
                if(headerMerger.hasReadGroupCollisions()) {
                    mappingToMerged.put(readGroup.getReadGroupId(),headerMerger.getReadGroupId(reader,readGroup.getReadGroupId()));
                    mergedToOriginalReadGroupMappings.put(headerMerger.getReadGroupId(reader,readGroup.getReadGroupId()),readGroup.getReadGroupId());
                } else {
                    mappingToMerged.put(readGroup.getReadGroupId(),readGroup.getReadGroupId());
                    mergedToOriginalReadGroupMappings.put(readGroup.getReadGroupId(),readGroup.getReadGroupId());
                }
            }

            originalToMergedReadGroupMappings.put(id,mappingToMerged);
        }

        if(enableLowMemorySharding) {
            for(SAMReaderID id: readerIDs) {
                File indexFile = findIndexFile(id.samFile);
                if(indexFile != null)
                    bamIndices.put(id,new GATKBAMIndex(indexFile));
            }
        }

        resourcePool.releaseReaders(readers);
    }

    /**
     * Returns Reads data structure containing information about the reads data sources placed in this pool as well as
     * information about how they are downsampled, sorted, and filtered
     * @return
     */
    public ReadProperties getReadsInfo() { return readProperties; }

    /**
     * Enable experimental low-memory sharding.
     * @param enable True to enable sharding.  False otherwise.
     */
    public void enableLowMemorySharding(final boolean enable) {
        enableLowMemorySharding = enable;
    }

    /**
     * Returns whether low-memory sharding is enabled.
     * @return True if enabled, false otherwise.
     */
    public boolean isLowMemoryShardingEnabled() {
        return enableLowMemorySharding;
    }

    /**
     * Checks to see whether any reads files are supplying data.
     * @return True if no reads files are supplying data to the traversal; false otherwise.
     */
    public boolean isEmpty() {
        return readProperties.getSAMReaderIDs().size() == 0;
    }

    /**
     * Gets the SAM file associated with a given reader ID.
     * @param id The reader for which to retrieve the source file.
     * @return the file actually associated with the id.
     */
    public File getSAMFile(SAMReaderID id) {
        return id.samFile;
    }

    /**
     * Returns readers used by this data source.
     * @return A list of SAM reader IDs.
     */
    public Collection<SAMReaderID> getReaderIDs() {
        return readerIDs;
    }

    /**
     * Retrieves the id of the reader which built the given read.
     * @param read The read to test.
     * @return ID of the reader.
     */
    public SAMReaderID getReaderID(SAMRecord read) {
        return resourcePool.getReaderID(read.getFileSource().getReader());
    }

    /**
     * Retrieves the current position within the BAM file.
     * @return A mapping of reader to current position.
     */
    public Map<SAMReaderID,SAMFileSpan> getCurrentPosition() {
        return readerPositions;
    }

    /**
     * Gets the merged header from the SAM file.
     * @return The merged header.
     */
    public SAMFileHeader getHeader() {
        return mergedHeader;
    }

    public SAMFileHeader getHeader(SAMReaderID id) {
        return resourcePool.getReadersWithoutLocking().getReader(id).getFileHeader();
    }

    /**
     * Gets the revised read group id mapped to this 'original' read group id.
     * @param reader for which to grab a read group.
     * @param originalReadGroupId ID of the original read group.
     * @return Merged read group ID.
     */
    public String getReadGroupId(final SAMReaderID reader, final String originalReadGroupId) {
        return originalToMergedReadGroupMappings.get(reader).get(originalReadGroupId);
    }

    /**
     * Gets the original read group id (as it was specified in the original input bam file) that maps onto
     * this 'merged' read group id.
     * @param mergedReadGroupId 'merged' ID of the read group (as it is presented by the read received from merged input stream).
     * @return Merged read group ID.
     */
    public String getOriginalReadGroupId(final String mergedReadGroupId) {
        return mergedToOriginalReadGroupMappings.get(mergedReadGroupId);
    }

    /**
     * No read group collisions at this time because only one SAM file is currently supported.
     * @return False always.
     */
    public boolean hasReadGroupCollisions() {
        return hasReadGroupCollisions;
    }

    /**
     * True if all readers have an index.
     * @return True if all readers have an index.
     */
    public boolean hasIndex() {
        if(enableLowMemorySharding)
            return readerIDs.size() == bamIndices.size();
        else {
            for(SAMFileReader reader: resourcePool.getReadersWithoutLocking()) {
                if(!reader.hasIndex())
                    return false;
            }
            return true;
        }
    }

    /**
     * Gets the index for a particular reader.  Always preloaded.
     * TODO: Should return object of type GATKBAMIndex, but cannot because there
     * TODO: is no parent class of both BAMIndex and GATKBAMIndex.  Change when new
     * TODO: sharding system goes live.      
     * @param id Id of the reader.
     * @return The index.  Will preload the index if necessary.
     */
    public Object getIndex(final SAMReaderID id) {
        if(enableLowMemorySharding)
            return bamIndices.get(id);
        else {
            SAMReaders readers = resourcePool.getReadersWithoutLocking();
            return readers.getReader(id).getBrowseableIndex();
        }
    }

    /**
     * Retrieves the sort order of the readers.
     * @return Sort order.  Can be unsorted, coordinate order, or query name order.
     */
    public SAMFileHeader.SortOrder getSortOrder() {
        return sortOrder;
    }

    /**
     * Gets the cumulative read metrics for shards already processed. 
     * @return Cumulative read metrics.
     */
    public ReadMetrics getCumulativeReadMetrics() {
        synchronized(readMetrics) {
            return readMetrics.clone();
        }
    }

    /**
     * Incorporate the given read metrics into the cumulative read metrics.
     * @param readMetrics The 'incremental' read metrics, to be incorporated into the cumulative metrics.
     */
    public void incorporateReadMetrics(final ReadMetrics readMetrics) {
        synchronized(this.readMetrics) {
            this.readMetrics.incrementMetrics(readMetrics);
        }
    }

    /**
     * Fill the given buffering shard with reads.
     * @param shard Shard to fill.
     * @return true if at the end of the stream.  False otherwise.
     */
    public void fillShard(Shard shard) {
        if(!shard.buffersReads())
            throw new ReviewedStingException("Attempting to fill a non-buffering shard.");

        SAMReaders readers = resourcePool.getAvailableReaders();
        // Cache the most recently viewed read so that we can check whether we've reached the end of a pair.
        SAMRecord read = null;

        CloseableIterator<SAMRecord> iterator = getIterator(readers,shard,sortOrder == SAMFileHeader.SortOrder.coordinate);
        while(!shard.isBufferFull() && iterator.hasNext()) {
            read = iterator.next();
            addReadToBufferingShard(shard,getReaderID(readers,read),read);
        }

        // If the reads are sorted in queryname order, ensure that all reads
        // having the same queryname become part of the same shard.
        if(sortOrder == SAMFileHeader.SortOrder.queryname) {
            while(iterator.hasNext()) {
                SAMRecord nextRead = iterator.next();
                if(read == null || !read.getReadName().equals(nextRead.getReadName()))
                    break;
                addReadToBufferingShard(shard,getReaderID(readers,nextRead),nextRead);
            }
        }

        iterator.close();
    }

    public StingSAMIterator seek(Shard shard) {
        // todo: refresh monolithic sharding implementation
        if(shard instanceof MonolithicShard)
            return seekMonolithic(shard);

        if(shard.buffersReads()) {
            return shard.iterator();
        }
        else {
            SAMReaders readers = resourcePool.getAvailableReaders();
            return getIterator(readers,shard,shard instanceof ReadShard);
        }
    }

    /**
     * Gets the reader associated with the given read.
     * @param readers Available readers.
     * @param read
     * @return
     */
    private SAMReaderID getReaderID(SAMReaders readers, SAMRecord read) {
        for(SAMReaderID id: getReaderIDs()) {
            if(readers.getReader(id) == read.getFileSource().getReader())
                return id;
        }
        throw new ReviewedStingException("Unable to find id for reader associated with read " + read.getReadName());
    }

    /**
     * Initialize the current reader positions
     * @param readers
     */
    private void initializeReaderPositions(SAMReaders readers) {
        for(SAMReaderID id: getReaderIDs())
            readerPositions.put(id,readers.getReader(id).getFilePointerSpanningReads());
    }

    /**
     * Get an iterator over the data types specified in the shard.
     * @param readers Readers from which to load data.
     * @param shard The shard specifying the data limits.
     * @param enableVerification True to verify.  For compatibility with old sharding strategy.
     *        TODO: Collapse this flag when the two sharding systems are merged.
     * @return An iterator over the selected data.
     */
    private StingSAMIterator getIterator(SAMReaders readers, Shard shard, boolean enableVerification) {
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate,readers.headers(),true);

        // Set up merging to dynamically merge together multiple BAMs.
        MergingSamRecordIterator mergingIterator = new MergingSamRecordIterator(headerMerger,readers.values(),true);

        for(SAMReaderID id: getReaderIDs()) {
            CloseableIterator<SAMRecord> iterator = null;
            if(!shard.isUnmapped() && shard.getFileSpans().get(id) == null)
                continue;
            iterator = shard.getFileSpans().get(id) != null ?
                    readers.getReader(id).iterator(shard.getFileSpans().get(id)) :
                    readers.getReader(id).queryUnmapped();
            if(readProperties.getReadBufferSize() != null)
                iterator = new BufferingReadIterator(iterator,readProperties.getReadBufferSize());
            if(shard.getGenomeLocs() != null)
                iterator = new IntervalOverlapFilteringIterator(iterator,shard.getGenomeLocs());
            mergingIterator.addIterator(readers.getReader(id),iterator);
        }

        return applyDecoratingIterators(shard.getReadMetrics(),
                enableVerification,
                readProperties.useOriginalBaseQualities(),
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(mergingIterator)),
                readProperties.getDownsamplingMethod().toFraction,
                readProperties.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                readProperties.getSupplementalFilters(),
                readProperties.getBAQCalculationMode(),
                readProperties.getBAQQualityMode(),
                readProperties.getRefReader(),
                readProperties.defaultBaseQualities());
    }

    /**
     * A stopgap measure to handle monolithic sharding
     * @param shard the (monolithic) shard.
     * @return An iterator over the monolithic shard.
     */
    private StingSAMIterator seekMonolithic(Shard shard) {
        SAMReaders readers = resourcePool.getAvailableReaders();

        // Set up merging and filtering to dynamically merge together multiple BAMs and filter out records not in the shard set.
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate,readers.headers(),true);
        MergingSamRecordIterator mergingIterator = new MergingSamRecordIterator(headerMerger,readers.values(),true);
        for(SAMReaderID id: getReaderIDs())
            mergingIterator.addIterator(readers.getReader(id),readers.getReader(id).iterator());

        return applyDecoratingIterators(shard.getReadMetrics(),
                shard instanceof ReadShard,
                readProperties.useOriginalBaseQualities(),
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(mergingIterator)),
                readProperties.getDownsamplingMethod().toFraction,
                readProperties.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                readProperties.getSupplementalFilters(),
                readProperties.getBAQCalculationMode(),
                readProperties.getBAQQualityMode(),
                readProperties.getRefReader(),
                readProperties.defaultBaseQualities());
    }

    /**
     * Adds this read to the given shard.
     * @param shard The shard to which to add the read.
     * @param id The id of the given reader.
     * @param read The read to add to the shard.
     */
    private void addReadToBufferingShard(Shard shard,SAMReaderID id,SAMRecord read) {
        SAMFileSpan endChunk = read.getFileSource().getFilePointer().getContentsFollowing();
        shard.addRead(read);
        readerPositions.put(id,endChunk);
    }

    /**
     * Filter reads based on user-specified criteria.
     *
     * @param readMetrics metrics to track when using this iterator.
     * @param enableVerification Verify the order of reads.
     * @param useOriginalBaseQualities True if original base qualities should be used.
     * @param wrappedIterator the raw data source.
     * @param downsamplingFraction whether and how much to downsample the reads themselves (not at a locus).
     * @param noValidationOfReadOrder Another trigger for the verifying iterator?  TODO: look into this.
     * @param supplementalFilters additional filters to apply to the reads.
     * @param defaultBaseQualities if the reads have incomplete quality scores, set them all to defaultBaseQuality.
     * @return An iterator wrapped with filters reflecting the passed-in parameters.  Will not be null.
     */
    protected StingSAMIterator applyDecoratingIterators(ReadMetrics readMetrics,
                                                        boolean enableVerification,
                                                        boolean useOriginalBaseQualities,
                                                        StingSAMIterator wrappedIterator,
                                                        Double downsamplingFraction,
                                                        Boolean noValidationOfReadOrder,
                                                        Collection<ReadFilter> supplementalFilters,
                                                        BAQ.CalculationMode cmode,
                                                        BAQ.QualityMode qmode,
                                                        IndexedFastaSequenceFile refReader,
                                                        byte defaultBaseQualities) {
        if ( useOriginalBaseQualities || defaultBaseQualities >= 0 )
            // only wrap if we are replacing the original qualitiies or using a default base quality
            wrappedIterator = new ReadFormattingIterator(wrappedIterator, useOriginalBaseQualities, defaultBaseQualities);

        // NOTE: this (and other filtering) should be done before on-the-fly sorting
        //  as there is no reason to sort something that we will end of throwing away
        if (downsamplingFraction != null)
            wrappedIterator = new DownsampleIterator(wrappedIterator, downsamplingFraction);

        // unless they've said not to validate read ordering (!noValidationOfReadOrder) and we've enabled verification,
        // verify the read ordering by applying a sort order iterator
        if (!noValidationOfReadOrder && enableVerification)
            wrappedIterator = new VerifyingSamIterator(genomeLocParser,wrappedIterator);

        if (cmode != BAQ.CalculationMode.OFF)
            wrappedIterator = new BAQSamIterator(refReader, wrappedIterator, cmode, qmode);

        wrappedIterator = StingSAMIteratorAdapter.adapt(new CountingFilteringIterator(readMetrics,wrappedIterator,supplementalFilters));

        return wrappedIterator;
    }

    private class SAMResourcePool {
        /**
         * How many entries can be cached in this resource pool?
         */
        private final int maxEntries;

        /**
         * All iterators of this reference-ordered data.
         */
        private List<SAMReaders> allResources = new ArrayList<SAMReaders>();

        /**
         * All iterators that are not currently in service.
         */
        private List<SAMReaders> availableResources = new ArrayList<SAMReaders>();

        public SAMResourcePool(final int maxEntries) {
            this.maxEntries = maxEntries;
        }

        /**
         * Dangerous internal method; retrieves any set of readers, whether in iteration or not.
         * Used to handle non-exclusive, stateless operations, such as index queries.
         * @return Any collection of SAMReaders, whether in iteration or not.
         */
        protected SAMReaders getReadersWithoutLocking() {
            synchronized(this) {
                if(allResources.size() == 0)
                    createNewResource();
            }
            return allResources.get(0);
        }

        /**
         * Choose a set of readers from the pool to use for this query.  When complete,
         * @return
         */
        public synchronized SAMReaders getAvailableReaders() {
            if(availableResources.size() == 0)
                createNewResource();
            SAMReaders readers = availableResources.get(0);
            availableResources.remove(readers);
            return readers;
        }

        public synchronized void releaseReaders(SAMReaders readers) {
            if(!allResources.contains(readers))
                throw new ReviewedStingException("Tried to return readers from the pool that didn't originate in the pool.");
            availableResources.add(readers);
        }

        /**
         * Gets the reader id for the given reader.
         * @param reader Reader for which to determine the id.
         * @return id of the given reader.
         */
        protected synchronized SAMReaderID getReaderID(SAMFileReader reader) {
            for(SAMReaders readers: allResources) {
                SAMReaderID id = readers.getReaderID(reader);
                if(id != null)
                    return id;
            }
            throw new ReviewedStingException("No such reader id is available");
        }

        private synchronized void createNewResource() {
            if(allResources.size() > maxEntries)
                throw new ReviewedStingException("Cannot create a new resource pool.  All resources are in use.");
            SAMReaders readers = new SAMReaders(readerIDs, validationStringency);
            allResources.add(readers);
            availableResources.add(readers);
        }

    }

    /**
     * A collection of readers derived from a reads metadata structure.
     */
    private class SAMReaders implements Iterable<SAMFileReader> {
        /**
         * Internal storage for a map of id -> reader.
         */
        private final Map<SAMReaderID,SAMFileReader> readers = new LinkedHashMap<SAMReaderID,SAMFileReader>();

        /**
         * Derive a new set of readers from the Reads metadata.
         * @param readerIDs reads to load.
         * @param validationStringency validation stringency.
         */
        public SAMReaders(Collection<SAMReaderID> readerIDs, SAMFileReader.ValidationStringency validationStringency) {
            for(SAMReaderID readerID: readerIDs) {
                SAMFileReader reader = new SAMFileReader(readerID.samFile);
                reader.setSAMRecordFactory(factory);
                reader.enableFileSource(true);
                reader.enableIndexMemoryMapping(false);
                if(!enableLowMemorySharding)
                    reader.enableIndexCaching(true);
                reader.setValidationStringency(validationStringency);

                final SAMFileHeader header = reader.getFileHeader();
                logger.debug(String.format("Sort order is: " + header.getSortOrder()));

                readers.put(readerID,reader);
            }
        }

        /**
         * Retrieve the reader from the data structure.
         * @param id The ID of the reader to retrieve.
         * @return the reader associated with the given id.
         */
        public SAMFileReader getReader(SAMReaderID id) {
            if(!readers.containsKey(id))
                throw new NoSuchElementException("No reader is associated with id " + id);
            return readers.get(id);
        }

        /**
         * Searches for the reader id of this reader.
         * @param reader Reader for which to search.
         * @return The id associated the given reader, or null if the reader is not present in this collection.
         */
        protected SAMReaderID getReaderID(SAMFileReader reader) {
            for(Map.Entry<SAMReaderID,SAMFileReader> entry: readers.entrySet()) {
                if(reader == entry.getValue())
                    return entry.getKey();
            }
            // Not found? return null.
            return null;
        }

        /**
         * Returns an iterator over all readers in this structure.
         * @return An iterator over readers.
         */
        public Iterator<SAMFileReader> iterator() {
            return readers.values().iterator();
        }

        /**
         * Returns whether any readers are present in this structure.
         * @return
         */
        public boolean isEmpty() {
            return readers.isEmpty();
        }

        /**
         * Gets all the actual readers out of this data structure.
         * @return A collection of the readers.
         */
        public Collection<SAMFileReader> values() {
            return readers.values();
        }

        /**
         * Gets all the actual readers out of this data structure.
         * @return A collection of the readers.
         */
        public Collection<SAMFileHeader> headers() {
            ArrayList<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(readers.size());
            for (SAMFileReader reader : values())
                headers.add(reader.getFileHeader());
            return headers;
        }
    }

    private class ReleasingIterator implements StingSAMIterator {
        /**
         * The resource acting as the source of the data.
         */
        private final SAMReaders resource;

        /**
         * The iterator to wrap.
         */
        private final StingSAMIterator wrappedIterator;

        public ReleasingIterator(SAMReaders resource, StingSAMIterator wrapped) {
            this.resource = resource;
            this.wrappedIterator = wrapped;
        }

        public ReleasingIterator iterator() {
            return this;
        }

        public void remove() {
            throw new UnsupportedOperationException("Can't remove from a StingSAMIterator");
        }

        public void close() {
            wrappedIterator.close();
            resourcePool.releaseReaders(resource);
        }

        public boolean hasNext() {
            return wrappedIterator.hasNext();
        }

        public SAMRecord next() {
            return wrappedIterator.next();
        }
    }

    /**
     * Maps read groups in the original SAMFileReaders to read groups in
     */
    private class ReadGroupMapping extends HashMap<String,String> {}

    /**
     * Filters out reads that do not overlap the current GenomeLoc.
     * Note the custom implementation: BAM index querying returns all reads that could
     * possibly overlap the given region (and quite a few extras).  In order not to drag
     * down performance, this implementation is highly customized to its task. 
     */
    private class IntervalOverlapFilteringIterator implements CloseableIterator<SAMRecord> {
        /**
         * The wrapped iterator.
         */
        private CloseableIterator<SAMRecord> iterator;

        /**
         * The next read, queued up and ready to go.
         */
        private SAMRecord nextRead;

        /**
         * Rather than using the straight genomic bounds, use filter out only mapped reads.
         */
        private boolean keepOnlyUnmappedReads;

        /**
         * Custom representation of interval bounds.
         * Makes it simpler to track current position. 
         */
        private int[] intervalContigIndices;
        private int[] intervalStarts;
        private int[] intervalEnds;

        /**
         * Position within the interval list.
         */
        private int currentBound = 0;

        public IntervalOverlapFilteringIterator(CloseableIterator<SAMRecord> iterator, List<GenomeLoc> intervals) {
            this.iterator = iterator;

            // Look at the interval list to detect whether we should worry about unmapped reads.
            // If we find a mix of mapped/unmapped intervals, throw an exception.
            boolean foundMappedIntervals = false;
            for(GenomeLoc location: intervals) {
                if(! GenomeLoc.isUnmapped(location))
                    foundMappedIntervals = true;
                keepOnlyUnmappedReads |= GenomeLoc.isUnmapped(location);
            }


            if(foundMappedIntervals) {
                if(keepOnlyUnmappedReads)
                    throw new ReviewedStingException("Tried to apply IntervalOverlapFilteringIterator to a mixed of mapped and unmapped intervals.  Please apply this filter to only mapped or only unmapped reads");
                this.intervalContigIndices = new int[intervals.size()];
                this.intervalStarts = new int[intervals.size()];
                this.intervalEnds = new int[intervals.size()];
                int i = 0;
                for(GenomeLoc interval: intervals) {
                    intervalContigIndices[i] = interval.getContigIndex();
                    intervalStarts[i] = interval.getStart();
                    intervalEnds[i] = interval.getStop();
                    i++;
                }
            }
            
            advance();
        }

        public boolean hasNext() {
            return nextRead != null;
        }

        public SAMRecord next() {
            if(nextRead == null)
                throw new NoSuchElementException("No more reads left in this iterator.");
            SAMRecord currentRead = nextRead;
            advance();
            return currentRead;
        }

        public void remove() {
            throw new UnsupportedOperationException("Cannot remove from an IntervalOverlapFilteringIterator");
        }


        public void close() {
            iterator.close();
        }

        private void advance() {
            nextRead = null;

            if(!iterator.hasNext())
                return;

            SAMRecord candidateRead = iterator.next();
            while(nextRead == null && (keepOnlyUnmappedReads || currentBound < intervalStarts.length)) {
                if(!keepOnlyUnmappedReads) {
                    // Mapped read filter; check against GenomeLoc-derived bounds.
                    if(readEndsOnOrAfterStartingBound(candidateRead)) {
                        // This read ends after the current interval begins.
                        // Promising, but this read must be checked against the ending bound.
                        if(readStartsOnOrBeforeEndingBound(candidateRead)) {
                            // Yes, this read is within both bounds.  This must be our next read.
                            nextRead = candidateRead;
                            break;
                        }
                        else {
                            // Oops, we're past the end bound.  Increment the current bound and try again.
                            currentBound++;
                            continue;
                        }
                    }
                }
                else {
                    // Found an unmapped read.  We're done.
                    if(candidateRead.getReadUnmappedFlag()) {
                        nextRead = candidateRead;
                        break;
                    }
                }

                // No more reads available.  Stop the search.
                if(!iterator.hasNext())
                    break;

                // No reasonable read found; advance the iterator.
                candidateRead = iterator.next();
            }
        }

        /**
         * Check whether the read lies after the start of the current bound.  If the read is unmapped but placed, its
         * end will be distorted, so rely only on the alignment start.
         * @param read The read to position-check.
         * @return True if the read starts after the current bounds.  False otherwise.
         */
        private boolean readEndsOnOrAfterStartingBound(final SAMRecord read) {
            return
                    // Read ends on a later contig, or...
                    read.getReferenceIndex() > intervalContigIndices[currentBound] ||
                    // Read ends of this contig...
                    (read.getReferenceIndex() == intervalContigIndices[currentBound] &&
                            // either after this location, or...
                            (read.getAlignmentEnd() >= intervalStarts[currentBound] ||
                            // read is unmapped but positioned and alignment start is on or after this start point.
                            (read.getReadUnmappedFlag() && read.getAlignmentStart() >= intervalStarts[currentBound])));
        }

        /**
         * Check whether the read lies before the end of the current bound.
         * @param read The read to position-check.
         * @return True if the read starts after the current bounds.  False otherwise.
         */
        private boolean readStartsOnOrBeforeEndingBound(final SAMRecord read) {
            return
                    // Read starts on a prior contig, or...
                    read.getReferenceIndex() < intervalContigIndices[currentBound] ||
                    // Read starts on this contig and the alignment start is registered before this end point.
                   (read.getReferenceIndex() == intervalContigIndices[currentBound] && read.getAlignmentStart() <= intervalEnds[currentBound]);
        }
    }

    /**
     * Locates the index file alongside the given BAM, if present.
     * TODO: This is currently a hachetjob that reaches into Picard and pulls out its index file locator.  Replace with something more permanent.
     * @param bamFile The data file to use.
     * @return A File object if the index file is present; null otherwise.
     */
    private File findIndexFile(File bamFile) {
        File indexFile;

        try {
            Class bamFileReaderClass = Class.forName("net.sf.samtools.BAMFileReader");
            Method indexFileLocator = bamFileReaderClass.getDeclaredMethod("findIndexFile",File.class);
            indexFileLocator.setAccessible(true);
            indexFile = (File)indexFileLocator.invoke(null,bamFile);
        }
        catch(ClassNotFoundException ex) {
            throw new ReviewedStingException("Unable to locate BAMFileReader class, used to check for index files");
        }
        catch(NoSuchMethodException ex) {
            throw new ReviewedStingException("Unable to locate Picard index file locator.");
        }
        catch(IllegalAccessException ex) {
            throw new ReviewedStingException("Unable to access Picard index file locator.");
        }
        catch(InvocationTargetException ex) {
            throw new ReviewedStingException("Unable to invoke Picard index file locator.");
        }

        return indexFile;
    }
}


