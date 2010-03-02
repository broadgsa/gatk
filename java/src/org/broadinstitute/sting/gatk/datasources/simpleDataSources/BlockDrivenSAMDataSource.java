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
    private final SAMFileHeader header;

    /**
     * Create a new block-aware SAM data source given the supplied read metadata.
     * @param reads The read metadata.
     */
    public BlockDrivenSAMDataSource(Reads reads) {
        super(reads);

        logger.warn("Experimental sharding is enabled.  Many use cases are not supported.  Please use with care.");

        resourcePool = new SAMResourcePool(Integer.MAX_VALUE);
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        header = new SamFileHeaderMerger(readers,SAMFileHeader.SortOrder.coordinate,true).getMergedHeader(); 
        resourcePool.releaseReaders(readers);
    }

    public boolean hasIndex() {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        try {
            return hasIndex(readers);
        }
        finally {
            resourcePool.releaseReaders(readers);
        }
    }

    /**
     * Report whether a given collection of SAM file readers is indexed.
     * @param readers The collection of readers.
     * @return True if the given collection of readers is indexed.
     */
    private boolean hasIndex(Collection<SAMFileReader> readers) {
        for(SAMFileReader reader: readers) {
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
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();

        try {
            if(readers.size() == 0)
                return Collections.emptyList();

            // All readers will have the same bin structure, so just use the first bin as an example.
            SAMFileReader2 reader = (SAMFileReader2)readers.iterator().next();
            return reader.getOverlappingBins(location.getContig(),(int)location.getStart(),(int)location.getStop());
        }
        finally {
            resourcePool.releaseReaders(readers);
        }
    }

    /**
     * Gets the file pointers bounded by this bin, grouped by the reader of origination.
     * @param bin The bin for which to load data.
     * @return A map of the file pointers bounding the bin.
     */
    public Map<SAMFileReader2,List<Chunk>> getFilePointersBounding(Bin bin) {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        try {
            Map<SAMFileReader2,List<Chunk>> filePointers = new HashMap<SAMFileReader2,List<Chunk>>();
            for(SAMFileReader reader: readers) {
                SAMFileReader2 reader2 = (SAMFileReader2)reader;
                filePointers.put(reader2,reader2.getFilePointersBounding(bin));
            }
            return filePointers;
        }
        finally {
            resourcePool.releaseReaders(readers);
        }
    }

    /**
     * Retrieves the current position within the BAM file.
     * @return A mapping of reader to current position.
     */
    public Map<SAMFileReader2,Chunk> getCurrentPosition() {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        try {
            Map<SAMFileReader2,Chunk> currentPositions = new HashMap<SAMFileReader2,Chunk>();
            for(SAMFileReader reader: readers) {
                SAMFileReader2 reader2 = (SAMFileReader2)reader;
                currentPositions.put(reader2,reader2.getCurrentPosition());
            }
            return currentPositions;
        }
        finally {
            resourcePool.releaseReaders(readers);
        }
    }

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    public int getNumIndexLevels() {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        try {
            if(readers.size() == 0)
                throw new StingException("Unable to determine number of index levels; no BAMs are present.");
            if(!hasIndex(readers))
                throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
            SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
            return firstReader.getNumIndexLevels();
        }
        finally {
            resourcePool.releaseReaders(readers);
        }            
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin) {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        try {
            if(readers.size() == 0)
                throw new StingException("Unable to determine number of level for bin; no BAMs are present.");
            if(!hasIndex(readers))
                throw new SAMException("Unable to determine number of level for bin; BAM file index is not present.");
            SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
            return firstReader.getLevelForBin(bin);
        }
        finally {
            resourcePool.releaseReaders(readers);
        }
    }

    /**
     * Gets the first locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getFirstLocusInBin(final Bin bin) {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        try {
            if(readers.size() == 0)
                throw new StingException("Unable to determine number of level for bin; no BAMs are present.");
            if(!hasIndex(readers))
                throw new SAMException("Unable to determine number of level for bin; BAM file index is not present.");
            SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
            return firstReader.getFirstLocusInBin(bin);
        }
        finally {
            resourcePool.releaseReaders(readers);
        }
    }

    /**
     * Gets the last locus that this bin can index into.
     * @param bin The bin to test.
     * @return The last position that the given bin can represent.
     */
    public int getLastLocusInBin(final Bin bin) {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();
        try {
            if(readers.size() == 0)
                throw new StingException("Unable to determine number of level for bin; no BAMs are present.");
            if(!hasIndex(readers))
                throw new SAMException("Unable to determine number of level for bin; BAM file index is not present.");
            SAMFileReader2 firstReader = (SAMFileReader2)readers.iterator().next();
            return firstReader.getLastLocusInBin(bin);
        }
        finally {
            resourcePool.releaseReaders(readers);
        }
    }

    /**
     * Fill the given buffering shard with reads.
     * @param shard Shard to fill.
     * @return true if at the end of the stream.  False otherwise.
     */
    public boolean fillShard(BAMFormatAwareShard shard) {
        if(!shard.buffersReads())
            throw new StingException("Attempting to fill a non-buffering shard.");

        // Since the beginning of time for the GATK, enableVerification has been true only for ReadShards.  I don't
        // know why this is.  Please add a comment here if you do.
        boolean enableVerification = shard instanceof ReadShard;

        CloseableIterator<SAMRecord> iterator = getIterator(shard,enableVerification);
        while(!shard.isBufferFull() && iterator.hasNext())
            shard.addRead(iterator.next());

        boolean atEndOfStream = !iterator.hasNext();

        iterator.close();

        return atEndOfStream;
    }

    public StingSAMIterator seek(Shard shard) {
        if(!(shard instanceof BAMFormatAwareShard))
            throw new StingException("BlockDrivenSAMDataSource cannot operate on shards of type: " + shard.getClass());
        BAMFormatAwareShard bamAwareShard = (BAMFormatAwareShard)shard;

        if(bamAwareShard.buffersReads()) {
            return bamAwareShard.iterator();    
        }
        else {
            // Since the beginning of time for the GATK, enableVerification has been true only for ReadShards.  I don't
            // know why this is.  Please add a comment here if you do.
            boolean enableVerification = shard instanceof ReadShard;
            return getIterator(bamAwareShard,enableVerification);
        }
    }

    private StingSAMIterator getIterator(BAMFormatAwareShard shard, boolean enableVerification) {
        Collection<SAMFileReader> readers = resourcePool.getAvailableReaders();        

        Map<SAMFileReader,CloseableIterator<SAMRecord>> readerToIteratorMap = new HashMap<SAMFileReader,CloseableIterator<SAMRecord>>();
        for(SAMFileReader reader: readers) {
            SAMFileReader2 reader2 = (SAMFileReader2)reader;
            List<Chunk> chunks = shard.getChunks().get(reader2);
            readerToIteratorMap.put(reader2,reader2.iterator(chunks));
        }

        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers,SAMFileHeader.SortOrder.coordinate,true);

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
     * Gets the merged header from the SAM file.
     * @return The merged header.
     */
    public SAMFileHeader getHeader() {
        return header;
    }

    /**
     * Currently unsupported.
     * @return
     */
    public Collection<SAMFileReader> getReaders() {
        throw new StingException("Currently unable to get readers for shard-based fields.");
    }

    /**
     * No read group collisions at this time because only one SAM file is currently supported.
     * @return False always.
     */
    public boolean hasReadGroupCollisions() {
        return false;
    }

    /**
     * Currently unsupported.
     * @return
     */
    public String getReadGroupId(final SAMFileReader reader, final String originalReadGroupId) {
        throw new UnsupportedOperationException("Getting read group ID from this experimental SAM reader is not currently supported.");
    }

    private class SAMResourcePool {
        /**
         * How many entries can be cached in this resource pool?
         */
        private final int maxEntries;

        /**
         * All iterators of this reference-ordered data.
         */
        private List<SAMFileReaders> allResources = new ArrayList<SAMFileReaders>();

        /**
         * All iterators that are not currently in service.
         */
        private List<SAMFileReaders> availableResources = new ArrayList<SAMFileReaders>();

        public SAMResourcePool(final int maxEntries) {
            this.maxEntries = maxEntries;
        }

        /**
         * Choose a set of readers from the pool to use for this query.  When complete,
         * @return
         */
        public synchronized Collection<SAMFileReader> getAvailableReaders() {
            if(availableResources.size() == 0)
                createNewResource();
            SAMFileReaders readers = availableResources.get(0);
            availableResources.remove(readers);
            return readers;
        }

        public synchronized void releaseReaders(Collection<SAMFileReader> readers) {
            if(!allResources.contains(readers))
                throw new StingException("Tried to return readers from the pool that didn't originate in the pool.");
            availableResources.add((SAMFileReaders)readers);
        }

        private synchronized void createNewResource() {
            if(allResources.size() > maxEntries)
                throw new StingException("Cannot create a new resource pool.  All resources are in use.");
            SAMFileReaders readers = new SAMFileReaders(reads);
            allResources.add(readers);
            availableResources.add(readers);
        }

        /**
         * A collection of readers derived from a reads metadata structure.
         */
        private class SAMFileReaders extends ArrayList<SAMFileReader> {
            /**
             * Derive a new set of readers from the Reads metadata.
             * @param sourceInfo Metadata for the reads to load.
             */
            public SAMFileReaders(Reads sourceInfo) {
                for(File readsFile: sourceInfo.getReadsFiles()) {
                    SAMFileReader2 reader = new SAMFileReader2(readsFile);
                    reader.setValidationStringency(sourceInfo.getValidationStringency());
                    add(reader);
                }
            }
        }        
    }

    private class ReleasingIterator implements StingSAMIterator {
        /**
         * The resource acting as the source of the data.
         */
        private final Collection<SAMFileReader> resource;

        /**
         * The iterator to wrap.
         */
        private final StingSAMIterator wrappedIterator;

        public Reads getSourceInfo() {
            return wrappedIterator.getSourceInfo();
        }

        public ReleasingIterator( Collection<SAMFileReader> resource, StingSAMIterator wrapped ) {
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
}
