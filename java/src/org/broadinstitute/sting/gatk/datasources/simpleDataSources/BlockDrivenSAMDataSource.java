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

    private final SamFileHeaderMerger headerMerger;

    /**
     * Create a new block-aware SAM data source given the supplied read metadata.
     * @param reads The read metadata.
     */
    public BlockDrivenSAMDataSource(Reads reads) {
        super(reads);

        logger.warn("Experimental sharding is enabled.  Many use cases are not supported.  Please use with care.");

        Collection<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        for(File readsFile: reads.getReadsFiles()) {
            SAMFileReader2 reader = new SAMFileReader2(readsFile);
            reader.setValidationStringency(reads.getValidationStringency());
            readers.add(reader);
        }

        this.headerMerger = new SamFileHeaderMerger(readers,SAMFileHeader.SortOrder.coordinate,true);
    }

    public boolean hasIndex() {
        for(SAMFileReader reader: headerMerger.getReaders()) {
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
        if(headerMerger.getReaders().size() == 0)
            return Collections.emptyList();

        // All readers will have the same bin structure, so just use the first bin as an example.
        SAMFileReader2 reader = (SAMFileReader2)headerMerger.getReaders().iterator().next();
        return reader.getOverlappingBins(location.getContig(),(int)location.getStart(),(int)location.getStop());
    }

    /**
     * Gets the file pointers bounded by this bin, grouped by the reader of origination.
     * @param bin The bin for which to load data.
     * @return A map of the file pointers bounding the bin.
     */
    public Map<SAMFileReader2,List<Chunk>> getFilePointersBounding(final Bin bin) {
        Map<SAMFileReader2,List<Chunk>> filePointers = new HashMap<SAMFileReader2,List<Chunk>>();
        for(SAMFileReader reader: headerMerger.getReaders()) {
            SAMFileReader2 reader2 = (SAMFileReader2)reader;
            filePointers.put(reader2,reader2.getFilePointersBounding(bin));
        }
        return filePointers;
    }


    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    public int getNumIndexLevels() {
        if(headerMerger.getReaders().size() == 0)
            throw new StingException("Unable to determine number of index levels; no BAMs are present.");
        if(!hasIndex())
            throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
        SAMFileReader2 firstReader = (SAMFileReader2)headerMerger.getReaders().iterator().next();
        return firstReader.getNumIndexLevels();
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin) {
        if(headerMerger.getReaders().size() == 0)
            throw new StingException("Unable to determine number of level for bin; no BAMs are present.");
        if(!hasIndex())
            throw new SAMException("Unable to determine number of level for bin; BAM file index is not present.");
        SAMFileReader2 firstReader = (SAMFileReader2)headerMerger.getReaders().iterator().next();
        return firstReader.getLevelForBin(bin);
    }

    public StingSAMIterator seek(Shard shard) {
        if(!(shard instanceof BAMFormatAwareShard))
            throw new StingException("BlockDrivenSAMDataSource cannot operate on shards of type: " + shard.getClass());
        BAMFormatAwareShard bamAwareShard = (BAMFormatAwareShard)shard;

        // Since the beginning of time for the GATK, enableVerification has been true only for ReadShards.  I don't
        // know why this is.  Please add a comment here if you do.
        boolean enableVerification = shard instanceof ReadShard;

        if(shard instanceof ReadShard && reads.getReadsFiles().size() > 1)
            throw new StingException("Experimental read sharding cannot handle multiple BAM files at this point.");

        Map<SAMFileReader,CloseableIterator<SAMRecord>> readerToIteratorMap = new HashMap<SAMFileReader,CloseableIterator<SAMRecord>>();
        for(Map.Entry<SAMFileReader2,List<Chunk>> chunksByReader: bamAwareShard.getChunks().entrySet()) {
            SAMFileReader2 reader = chunksByReader.getKey();
            List<Chunk> chunks = chunksByReader.getValue();
            readerToIteratorMap.put(reader,reader.iterator(chunks));
        }

        // Set up merging and filtering to dynamically merge together multiple BAMs and filter out records not in the shard set.
        MergingSamRecordIterator mergingIterator = new MergingSamRecordIterator(headerMerger,readerToIteratorMap,true);
        FilteringIterator filteringIterator = new FilteringIterator(mergingIterator,new IntervalOverlappingFilter(shard.getGenomeLocs()));

        return applyDecoratingIterators(enableVerification,
                StingSAMIteratorAdapter.adapt(reads,filteringIterator),
                reads.getDownsamplingFraction(),
                reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                reads.getSupplementalFilters());
    }

    /**
     * Gets the merged header from the SAM file.
     * @return The merged header.
     */
    public SAMFileHeader getHeader() {
        return headerMerger.getMergedHeader();
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
}
