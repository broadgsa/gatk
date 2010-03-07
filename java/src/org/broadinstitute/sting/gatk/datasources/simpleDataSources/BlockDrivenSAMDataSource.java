package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.datasources.shards.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.filter.FilteringIterator;

import java.util.*;
import java.io.File;

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
    private final Map<SAMReaderID,Chunk> readerPositions = new HashMap<SAMReaderID,Chunk>();

    /**
     * Create a new block-aware SAM data source given the supplied read metadata.
     * @param reads The read metadata.
     */
    public BlockDrivenSAMDataSource(Reads reads) {
        super(reads);

        logger.warn("Experimental sharding is enabled.  Many use cases are not supported.  Please use with care.");

        resourcePool = new SAMResourcePool(Integer.MAX_VALUE);
        SAMReaders readers = resourcePool.getAvailableReaders();

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
     * Gets a list of the bins in each BAM file that overlap with the given interval list.
     * @param location Location for which to determine the bin.
     * @return A map of reader back to bin.
     */
    public List<Bin> getOverlappingBins(final GenomeLoc location) {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        if(readers.isEmpty())
            return Collections.emptyList();

        // All readers will have the same bin structure, so just use the first bin as an example.
        SAMFileReader2 reader = (SAMFileReader2)readers.iterator().next();
        return reader.getOverlappingBins(location.getContig(),(int)location.getStart(),(int)location.getStop());
    }

    /**
     * Gets the file pointers bounded by this bin, grouped by the reader of origination.
     * @param bin The bin for which to load data.
     * @return A map of the file pointers bounding the bin.
     */
    public Map<SAMReaderID,List<Chunk>> getFilePointersBounding(Bin bin) {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        Map<SAMReaderID,List<Chunk>> filePointers = new HashMap<SAMReaderID,List<Chunk>>();
        for(SAMReaderID id: getReaderIDs()) {
            SAMFileReader2 reader2 = (SAMFileReader2)readers.getReader(id);
            if(bin != null)
                filePointers.put(id,reader2.getFilePointersBounding(bin));
            else
                filePointers.put(id,Collections.<Chunk>emptyList());
        }
        return filePointers;
    }

    /**
     * Retrieves the current position within the BAM file.
     * @return A mapping of reader to current position.
     */
    public Map<SAMReaderID,Chunk> getCurrentPosition() {
        return readerPositions;
    }

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    public int getNumIndexLevels() {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        if(readers.isEmpty())
            throw new StingException("Unable to determine number of index levels; no BAMs are present.");
        if(!hasIndex())
            throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
        SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
        return firstReader.getNumIndexLevels();
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin) {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        if(readers.isEmpty())
            throw new StingException("Unable to determine number of level for bin; no BAMs are present.");
        if(!hasIndex())
            throw new SAMException("Unable to determine number of level for bin; BAM file index is not present.");
        SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
        return firstReader.getLevelForBin(bin);
    }

    /**
     * Gets the first locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getFirstLocusInBin(final Bin bin) {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        if(readers.isEmpty())
            throw new StingException("Unable to determine number of level for bin; no BAMs are present.");
        if(!hasIndex())
            throw new SAMException("Unable to determine number of level for bin; BAM file index is not present.");
        SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
        return firstReader.getFirstLocusInBin(bin);
    }

    /**
     * Gets the last locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getLastLocusInBin(final Bin bin) {
        SAMReaders readers = resourcePool.getReadersWithoutLocking();
        if(readers.isEmpty())
            throw new StingException("Unable to determine number of level for bin; no BAMs are present.");
        if(!hasIndex())
            throw new SAMException("Unable to determine number of level for bin; BAM file index is not present.");
        SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
        return firstReader.getLastLocusInBin(bin);
    }

    /**
     * Fill the given buffering shard with reads.
     * @param shard Shard to fill.
     * @return true if at the end of the stream.  False otherwise.
     */
    public void fillShard(BAMFormatAwareShard shard) {
        if(!shard.buffersReads())
            throw new StingException("Attempting to fill a non-buffering shard.");

        // Since the beginning of time for the GATK, enableVerification has been true only for ReadShards.  I don't
        // know why this is.  Please add a comment here if you do.
        boolean enableVerification = shard instanceof ReadShard;

        SAMReaders readers = resourcePool.getAvailableReaders();

        CloseableIterator<SAMRecord> iterator = getIterator(readers,shard,enableVerification);
        while(!shard.isBufferFull() && iterator.hasNext()) {
            SAMRecord read = iterator.next();
            Chunk endChunk = new Chunk(read.getCoordinates().getChunkEnd(),Long.MAX_VALUE);
            shard.addRead(read);
            readerPositions.put(getReaderID(readers,read),endChunk);
        }

        iterator.close();
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

    private void initializeReaderPositions(SAMReaders readers) {
        for(SAMReaderID id: getReaderIDs()) {
            SAMFileReader2 reader2 = (SAMFileReader2)readers.getReader(id);
            readerPositions.put(id,reader2.getCurrentPosition());
        }
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

    private StingSAMIterator getIterator(SAMReaders readers, BAMFormatAwareShard shard, boolean enableVerification) {
        Map<SAMFileReader,CloseableIterator<SAMRecord>> readerToIteratorMap = new HashMap<SAMFileReader,CloseableIterator<SAMRecord>>();
        for(SAMReaderID id: getReaderIDs()) {
            SAMFileReader2 reader2 = (SAMFileReader2)readers.getReader(id);
            List<Chunk> chunks = shard.getChunks().get(id);
            readerToIteratorMap.put(reader2,reader2.iterator(chunks));
        }

        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers.values(),SAMFileHeader.SortOrder.coordinate,true);

        // Set up merging and filtering to dynamically merge together multiple BAMs and filter out records not in the shard set.
        CloseableIterator<SAMRecord> iterator = new MergingSamRecordIterator(headerMerger,readerToIteratorMap,true);
        if(shard.getFilter() != null)
            iterator = new FilteringIterator(iterator,shard.getFilter());

        return applyDecoratingIterators(enableVerification,
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(reads,iterator)),
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

        Map<SAMFileReader,CloseableIterator<SAMRecord>> readerToIteratorMap = new HashMap<SAMFileReader,CloseableIterator<SAMRecord>>();
        for(SAMReaderID id: getReaderIDs()) {
            SAMFileReader2 reader2 = (SAMFileReader2)readers.getReader(id);
            readerToIteratorMap.put(reader2,reader2.iterator());
        }

        // Set up merging and filtering to dynamically merge together multiple BAMs and filter out records not in the shard set.
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers.values(),SAMFileHeader.SortOrder.coordinate,true);
        CloseableIterator<SAMRecord> iterator = new MergingSamRecordIterator(headerMerger,readerToIteratorMap,true);

        return applyDecoratingIterators(shard instanceof ReadShard,
                new ReleasingIterator(readers,StingSAMIteratorAdapter.adapt(reads,iterator)),
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
                SAMFileReader2 reader = new SAMFileReader2(readsFile,true);
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
