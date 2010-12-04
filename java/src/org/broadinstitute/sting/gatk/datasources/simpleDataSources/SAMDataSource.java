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

package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.DownsamplingMethod;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.BAMFormatAwareShard;
import org.broadinstitute.sting.gatk.datasources.shards.MonolithicShard;
import org.broadinstitute.sting.gatk.datasources.shards.ReadShard;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.filters.CountingFilteringIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.BAQ;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.util.*;

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 2:36:16 PM
 * <p/>
 * Converts shards to SAM iterators over the specified region
 */
public class SAMDataSource implements SimpleDataSource {
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
    private final List<SAMReaderID> readerIDs;

    /**
     * How strict are the readers driving this data source.
     */
    private final SAMFileReader.ValidationStringency validationStringency;

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
     * Create a new SAM data source given the supplied read metadata.
     * @param samFiles list of reads files.
     */
    public SAMDataSource(List<SAMReaderID> samFiles,GenomeLocParser genomeLocParser) {
        this(
                samFiles,
                genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<SamRecordFilter>(),
                false,
                false
        );
    }

    /**
     * See complete constructor.  Does not enable BAQ by default.
     */
    public SAMDataSource(
            List<SAMReaderID> samFiles,
            GenomeLocParser genomeLocParser,
            boolean useOriginalBaseQualities,
            SAMFileReader.ValidationStringency strictness,
            Integer readBufferSize,
            DownsamplingMethod downsamplingMethod,
            ValidationExclusion exclusionList,
            Collection<SamRecordFilter> supplementalFilters,
            boolean includeReadsWithDeletionAtLoci,
            boolean generateExtendedEvents ) {
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
                BAQ.Mode.NONE, null                 // no BAQ
                );
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
     */
    public SAMDataSource(
            List<SAMReaderID> samFiles,
            GenomeLocParser genomeLocParser,
            boolean useOriginalBaseQualities,
            SAMFileReader.ValidationStringency strictness,
            Integer readBufferSize,
            DownsamplingMethod downsamplingMethod,
            ValidationExclusion exclusionList,
            Collection<SamRecordFilter> supplementalFilters,
            boolean includeReadsWithDeletionAtLoci,
            boolean generateExtendedEvents,
            BAQ.Mode Mode,
            IndexedFastaSequenceFile refReader
    ) {
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
                Mode, refReader );
        
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

        resourcePool.releaseReaders(readers);
    }

    /**
     * Returns Reads data structure containing information about the reads data sources placed in this pool as well as
     * information about how they are downsampled, sorted, and filtered
     * @return
     */
    public ReadProperties getReadsInfo() { return readProperties; }

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
    public List<SAMReaderID> getReaderIDs() {
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
     * @return
     */
    public boolean hasIndex() {
        for(SAMFileReader reader: resourcePool.getReadersWithoutLocking()) {
            if(!reader.hasIndex())
                return false;
        }
        return true;
    }

    /**
     * Gets the index for a particular reader.  Always preloaded.
     * @param id Id of the reader.
     * @return The index.  Will preload the index if necessary.
     */
    public BrowseableBAMIndex getIndex(final SAMReaderID id) {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        return readers.getReader(id).getBrowseableIndex();
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
    public void fillShard(BAMFormatAwareShard shard) {
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

        if(!(shard instanceof BAMFormatAwareShard))
            throw new ReviewedStingException("BlockDrivenSAMDataSource cannot operate on shards of type: " + shard.getClass());
        BAMFormatAwareShard bamAwareShard = (BAMFormatAwareShard)shard;

        if(bamAwareShard.buffersReads()) {
            return bamAwareShard.iterator();
        }
        else {
            SAMReaders readers = resourcePool.getAvailableReaders();
            return getIterator(readers,bamAwareShard,shard instanceof ReadShard);
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
    private StingSAMIterator getIterator(SAMReaders readers, BAMFormatAwareShard shard, boolean enableVerification) {
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate,readers.headers(),true);

        // Set up merging to dynamically merge together multiple BAMs.
        MergingSamRecordIterator mergingIterator = new MergingSamRecordIterator(headerMerger,readers.values(),true);

        for(SAMReaderID id: getReaderIDs()) {
            if(shard.getFileSpans().get(id) == null)
                continue;
            CloseableIterator<SAMRecord> iterator = readers.getReader(id).iterator(shard.getFileSpans().get(id));            
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
                readProperties.getBAQMode(), readProperties.getRefReader());
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
                readProperties.getBAQMode(), readProperties.getRefReader());
    }

    /**
     * Adds this read to the given shard.
     * @param shard The shard to which to add the read.
     * @param id The id of the given reader.
     * @param read The read to add to the shard.
     */
    private void addReadToBufferingShard(BAMFormatAwareShard shard,SAMReaderID id,SAMRecord read) {
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
     * @return An iterator wrapped with filters reflecting the passed-in parameters.  Will not be null.
     */
    protected StingSAMIterator applyDecoratingIterators(ReadMetrics readMetrics,
                                                        boolean enableVerification,
                                                        boolean useOriginalBaseQualities,
                                                        StingSAMIterator wrappedIterator,
                                                        Double downsamplingFraction,
                                                        Boolean noValidationOfReadOrder,
                                                        Collection<SamRecordFilter> supplementalFilters,
                                                        BAQ.Mode mode, IndexedFastaSequenceFile refReader ) {
        wrappedIterator = new ReadFormattingIterator(wrappedIterator, useOriginalBaseQualities);

        // NOTE: this (and other filtering) should be done before on-the-fly sorting
        //  as there is no reason to sort something that we will end of throwing away
        if (downsamplingFraction != null)
            wrappedIterator = new DownsampleIterator(wrappedIterator, downsamplingFraction);

        // unless they've said not to validate read ordering (!noValidationOfReadOrder) and we've enabled verification,
        // verify the read ordering by applying a sort order iterator
        if (!noValidationOfReadOrder && enableVerification)
            wrappedIterator = new VerifyingSamIterator(genomeLocParser,wrappedIterator);

        if (mode != BAQ.Mode.NONE)
            wrappedIterator = new BAQSamIterator(refReader, wrappedIterator, mode);

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
                reader.enableFileSource(true);
                reader.enableIndexCaching(true);
                reader.setValidationStringency(validationStringency);

                // If no read group is present, hallucinate one.
                // TODO: Straw poll to see whether this is really required.
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
         * Custom representation of interval bounds.
         * Makes it simpler to track current position. 
         */
        private int[] intervalStarts;
        private int[] intervalEnds;

        /**
         * Position within the interval list.
         */
        private int currentBound = 0;

        public IntervalOverlapFilteringIterator(CloseableIterator<SAMRecord> iterator, List<GenomeLoc> intervals) {
            this.iterator = iterator;
            this.intervalStarts = new int[intervals.size()];
            this.intervalEnds = new int[intervals.size()];
            int i = 0;
            for(GenomeLoc interval: intervals) {
                intervalStarts[i] = (int)interval.getStart();
                intervalEnds[i] = (int)interval.getStop();
                i++;
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
            while(nextRead == null && currentBound < intervalStarts.length) {
                if(candidateRead.getAlignmentEnd() >= intervalStarts[currentBound] ||
                  (candidateRead.getReadUnmappedFlag() && candidateRead.getAlignmentStart() >= intervalStarts[currentBound])) {
                    // This read ends after the current interval begins (or, if unmapped, starts within the bounds of the interval.
                    // Promising, but this read must be checked against the ending bound.
                    if(candidateRead.getAlignmentStart() <= intervalEnds[currentBound]) {
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

                // No more reads available.  Stop the search.
                if(!iterator.hasNext())
                    break;

                // No reasonable read found; advance the iterator.
                candidateRead = iterator.next();
            }
        }
    }
}


