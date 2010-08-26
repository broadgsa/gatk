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
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.picard.sam.MergingSamRecordIterator;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.BAMFormatAwareShard;
import org.broadinstitute.sting.gatk.datasources.shards.MonolithicShard;
import org.broadinstitute.sting.gatk.datasources.shards.ReadShard;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.filters.CountingFilteringIterator;
import org.broadinstitute.sting.utils.StingException;

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
    protected final ReadMetrics readMetrics;

    /**
     * Identifiers for the readers driving this data source.
     */
    protected final List<SAMReaderID> readerIDs;

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
     * Maps the SAM readers' original read group ids to their revised ids.
     */
    private final Map<SAMReaderID,ReadGroupMapping> mergedReadGroupMappings = new HashMap<SAMReaderID,ReadGroupMapping>();

    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(SAMDataSource.class);

    /**
     * A collection of readers driving the merging process.
     */
    private final SAMResourcePool resourcePool;

    /**
     * Create a new SAM data source given the supplied read metadata.
     * @param reads The read metadata.
     */
    public SAMDataSource(ReadProperties reads) {
        this.readProperties = reads;
        this.readMetrics = new ReadMetrics();

        readerIDs = reads.getSAMReaderIDs();
        for (SAMReaderID readerID : reads.getSAMReaderIDs()) {
            if (!readerID.samFile.canRead())
                throw new SimpleDataSourceLoadException("SAMDataSource: Unable to load file: " + readerID.samFile.getName());
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
                throw new StingException(String.format("Attempted to process mixed of files sorted as %s and %s.",this.sortOrder,sortOrder));

            // Update the sort order.
            this.sortOrder = sortOrder;
        }

        initializeReaderPositions(readers);

        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers.values(),SAMFileHeader.SortOrder.coordinate,true);
        mergedHeader = headerMerger.getMergedHeader();
        hasReadGroupCollisions = headerMerger.hasReadGroupCollisions();

        // cache the read group id (original) -> read group id (merged) mapping.
        for(SAMReaderID id: readerIDs) {
            SAMFileReader reader = readers.getReader(id);
            ReadGroupMapping mapping = new ReadGroupMapping();

            List<SAMReadGroupRecord> readGroups = reader.getFileHeader().getReadGroups();
            for(SAMReadGroupRecord readGroup: readGroups) {
                if(headerMerger.hasReadGroupCollisions())
                    mapping.put(readGroup.getReadGroupId(),headerMerger.getReadGroupId(reader,readGroup.getReadGroupId()));
                else
                    mapping.put(readGroup.getReadGroupId(),readGroup.getReadGroupId());
            }

            mergedReadGroupMappings.put(id,mapping);
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
        return mergedReadGroupMappings.get(reader).get(originalReadGroupId);
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
            throw new StingException("Attempting to fill a non-buffering shard.");

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
            throw new StingException("BlockDrivenSAMDataSource cannot operate on shards of type: " + shard.getClass());
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
        throw new StingException("Unable to find id for reader associated with read " + read.getReadName());
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
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers.values(),SAMFileHeader.SortOrder.coordinate,true);

        // Set up merging to dynamically merge together multiple BAMs.
        MergingSamRecordIterator mergingIterator = new MergingSamRecordIterator(headerMerger,true);
        for(SAMReaderID id: getReaderIDs()) {
            if(shard.getFileSpans().get(id) == null)
                continue;
            CloseableIterator<SAMRecord> iterator = readers.getReader(id).iterator(shard.getFileSpans().get(id));
            if(readProperties.getReadBufferSize() != null)
                iterator = new BufferingReadIterator(iterator,readProperties.getReadBufferSize());
            if(shard.getFilter() != null)
                iterator = new FilteringIterator(iterator,shard.getFilter()); // not a counting iterator because we don't want to show the filtering of reads
            mergingIterator.addIterator(readers.getReader(id),iterator);
        }

        return applyDecoratingIterators(shard.getReadMetrics(),
                enableVerification,
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(mergingIterator)),
                readProperties.getDownsamplingMethod().toFraction,
                readProperties.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                readProperties.getSupplementalFilters());
    }

    /**
     * A stopgap measure to handle monolithic sharding
     * @param shard the (monolithic) shard.
     * @return An iterator over the monolithic shard.
     */
    private StingSAMIterator seekMonolithic(Shard shard) {
        SAMReaders readers = resourcePool.getAvailableReaders();

        // Set up merging and filtering to dynamically merge together multiple BAMs and filter out records not in the shard set.
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers.values(),SAMFileHeader.SortOrder.coordinate,true);
        MergingSamRecordIterator mergingIterator = new MergingSamRecordIterator(headerMerger,true);
        for(SAMReaderID id: getReaderIDs())
            mergingIterator.addIterator(readers.getReader(id),readers.getReader(id).iterator());

        return applyDecoratingIterators(shard.getReadMetrics(),
                shard instanceof ReadShard,
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(mergingIterator)),
                readProperties.getDownsamplingMethod().toFraction,
                readProperties.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                readProperties.getSupplementalFilters());
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
     * @param wrappedIterator the raw data source.
     * @param downsamplingFraction whether and how much to downsample the reads themselves (not at a locus).
     * @param noValidationOfReadOrder Another trigger for the verifying iterator?  TODO: look into this.
     * @param supplementalFilters additional filters to apply to the reads.
     * @return An iterator wrapped with filters reflecting the passed-in parameters.  Will not be null.
     */
    protected StingSAMIterator applyDecoratingIterators(ReadMetrics readMetrics,
                                                        boolean enableVerification,
                                                        StingSAMIterator wrappedIterator,
                                                        Double downsamplingFraction,
                                                        Boolean noValidationOfReadOrder,
                                                        Collection<SamRecordFilter> supplementalFilters) {
        wrappedIterator = new ReadFormattingIterator(wrappedIterator);

        // NOTE: this (and other filtering) should be done before on-the-fly sorting
        //  as there is no reason to sort something that we will end of throwing away
        if (downsamplingFraction != null)
            wrappedIterator = new DownsampleIterator(wrappedIterator, downsamplingFraction);

        // unless they've said not to validate read ordering (!noValidationOfReadOrder) and we've enabled verification,
        // verify the read ordering by applying a sort order iterator
        if (!noValidationOfReadOrder && enableVerification)
            wrappedIterator = new VerifyingSamIterator(wrappedIterator);

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
                throw new StingException("Tried to return readers from the pool that didn't originate in the pool.");
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
            throw new StingException("No such reader id is available");
        }

        private synchronized void createNewResource() {
            if(allResources.size() > maxEntries)
                throw new StingException("Cannot create a new resource pool.  All resources are in use.");
            SAMReaders readers = new SAMReaders(readProperties);
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
         * @param sourceInfo Metadata for the reads to load.
         */
        public SAMReaders(ReadProperties sourceInfo) {
            for(SAMReaderID readerID: sourceInfo.getSAMReaderIDs()) {
                SAMFileReader reader = new SAMFileReader(readerID.samFile);
                reader.enableFileSource(true);
                reader.enableIndexCaching(true);
                reader.setValidationStringency(sourceInfo.getValidationStringency());

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
}


