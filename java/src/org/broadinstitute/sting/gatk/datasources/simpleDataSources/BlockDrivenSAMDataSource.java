package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.datasources.shards.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.utils.StingException;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.filter.FilteringIterator;

import java.util.*;
import java.io.File;

import static net.sf.samtools.SAMFileHeader.SortOrder;

/**
 * An iterator that's aware of how data is stored on disk in SAM format.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockDrivenSAMDataSource extends SAMDataSource {
    /**
     * A collection of readers driving the merging process.
     */
    private final SAMResourcePool resourcePool;

    /**
     * The sort order of the BAM files.  Files without a sort order tag are assumed to be
     * in coordinate order.
     */
    private SAMFileHeader.SortOrder sortOrder = null;

    /**
     * The merged header.
     */
    private final SAMFileHeader mergedHeader;

    /**
     * Whether the read groups in overlapping files collide.
     */
    private final boolean hasReadGroupCollisions;

    /**
     * Maps the SAM readers' original read group ids to their revised ids.
     */
    private final Map<SAMReaderID,ReadGroupMapping> mergedReadGroupMappings = new HashMap<SAMReaderID,ReadGroupMapping>();

    /**
     * How far along is each reader?
     */
    private final Map<SAMReaderID,SAMFileSpan> readerPositions = new HashMap<SAMReaderID,SAMFileSpan>();

    /**
     * Create a new block-aware SAM data source given the supplied read metadata.
     * @param reads The read metadata.
     */
    public BlockDrivenSAMDataSource(Reads reads) {
        super(reads);

        resourcePool = new SAMResourcePool(Integer.MAX_VALUE);
        SAMReaders readers = resourcePool.getAvailableReaders();

        // Determine the sort order.
        for(SAMFileReader reader: readers.values()) {
            // Get the sort order, forcing it to coordinate if unsorted.
            SAMFileHeader header = reader.getFileHeader();
            SortOrder sortOrder = header.getSortOrder() != SortOrder.unsorted ? header.getSortOrder() : SortOrder.coordinate;

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
     * Retrieves the sort order of the readers.
     * @return Sort order.  Can be unsorted, coordinate order, or query name order.
     */
    public SortOrder getSortOrder() {
        return sortOrder;    
    }

    /**
     * Gets the index for a particular reader.  Always preloaded.
     * @param id Id of the reader.
     * @return The index.  Will preload the index if necessary.
     */
    public CachingBAMFileIndex getIndex(final SAMReaderID id) {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        return readers.getReader(id).getIndex(CachingBAMFileIndex.class);
    }

    /**
     * Retrieves the current position within the BAM file.
     * @return A mapping of reader to current position.
     */
    public Map<SAMReaderID,SAMFileSpan> getCurrentPosition() {
        return readerPositions;
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

        CloseableIterator<SAMRecord> iterator = getIterator(readers,shard,sortOrder == SortOrder.coordinate);
        while(!shard.isBufferFull() && iterator.hasNext()) {
            read = iterator.next();
            addReadToBufferingShard(shard,getReaderID(readers,read),read);
        }

        // If the reads are sorted in queryname order, ensure that all reads
        // having the same queryname become part of the same shard.
        if(sortOrder == SortOrder.queryname) {
            while(iterator.hasNext()) {
                SAMRecord nextRead = iterator.next();
                if(read == null || !read.getReadName().equals(nextRead.getReadName()))
                    break;    
                addReadToBufferingShard(shard,getReaderID(readers,nextRead),nextRead);
            }
        }

        iterator.close();
    }

    /**
     * Retrieves the id of the reader which built the given read.
     * @param read The read to test.
     * @return ID of the reader.
     */
    public SAMReaderID getReaderID(SAMRecord read) {
        return resourcePool.getReaderID(read.getReader());
    }    

    /**
     * Adds this read to the given shard.
     * @param shard The shard to which to add the read.
     * @param id The id of the given reader.
     * @param read The read to add to the shard.
     */
    private void addReadToBufferingShard(BAMFormatAwareShard shard,SAMReaderID id,SAMRecord read) {
        SAMFileSpan endChunk = read.getFilePointer().getContentsFollowing();
        shard.addRead(read);
        readerPositions.put(id,endChunk);
    }

    /**
     * Gets the reader associated with the given read.
     * @param readers Available readers.
     * @param read
     * @return
     */
    private SAMReaderID getReaderID(SAMReaders readers, SAMRecord read) {
        for(SAMReaderID id: getReaderIDs()) {
            if(readers.getReader(id) == read.getReader())
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
            if(shard.getFilter() != null)
                iterator = new FilteringIterator(iterator,shard.getFilter());
            mergingIterator.addIterator(readers.getReader(id),iterator);
        }

        return applyDecoratingIterators(enableVerification,
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(reads,mergingIterator)),
                reads.getDownsamplingFraction(),
                reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                reads.getSupplementalFilters());        
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

        return applyDecoratingIterators(shard instanceof ReadShard,
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(reads,mergingIterator)),
                reads.getDownsamplingFraction(),
                reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                reads.getSupplementalFilters());
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
     * No read group collisions at this time because only one SAM file is currently supported.
     * @return False always.
     */
    public boolean hasReadGroupCollisions() {
        return hasReadGroupCollisions;
    }

    /**
     * Gets the revised read group id mapped to this 'original' read group id.
     * @return Merged read group ID.
     */
    public String getReadGroupId(final SAMReaderID reader, final String originalReadGroupId) {
        return mergedReadGroupMappings.get(reader).get(originalReadGroupId);
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
            SAMReaders readers = new SAMReaders(reads);
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
        public SAMReaders(Reads sourceInfo) {
            for(File readsFile: sourceInfo.getReadsFiles()) {
                SAMFileReader reader = new SAMFileReader(readsFile,CachingBAMFileIndex.class,true);
                reader.setValidationStringency(sourceInfo.getValidationStringency());

                // If no read group is present, hallucinate one.
                // TODO: Straw poll to see whether this is really required.
                final SAMFileHeader header = reader.getFileHeader();
                logger.debug(String.format("Sort order is: " + header.getSortOrder()));

                if (reader.getFileHeader().getReadGroups().size() < 1) {
                    SAMReadGroupRecord rec = new SAMReadGroupRecord(readsFile.getName());
                    rec.setLibrary(readsFile.getName());
                    rec.setSample(readsFile.getName());

                    reader.getFileHeader().addReadGroup(rec);
                }

                readers.put(new SAMReaderID(readsFile),reader);
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

        public Reads getSourceInfo() {
            return wrappedIterator.getSourceInfo();
        }

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
