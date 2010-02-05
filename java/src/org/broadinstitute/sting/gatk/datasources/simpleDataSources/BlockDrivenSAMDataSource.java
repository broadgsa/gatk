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

import java.util.Collection;
import java.util.List;
import java.io.File;

/**
 * An iterator that's aware of how data is stored on disk in SAM format.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockDrivenSAMDataSource extends SAMDataSource {

    private final SAMFileReader2 reader;
    /**
     * Create a new block-aware SAM data source given the supplied read metadata.
     * @param reads The read metadata.
     */
    public BlockDrivenSAMDataSource(Reads reads) {
        super(reads);

        logger.warn("Experimental sharding is enabled.  Many use cases are not supported.  Please use with care.");

        if(reads.getReadsFiles().size() > 1)
            throw new StingException("Experimental sharding strategy cannot handle multiple BAM files at this point.");

        File readsFile = reads.getReadsFiles().get(0);
        reader = new SAMFileReader2(readsFile);
        reader.setValidationStringency(reads.getValidationStringency());
    }

    public boolean hasIndex() {
        return reader.hasIndex();
    }

    public List<Bin> getOverlappingBins(GenomeLoc location) {
        return reader.getOverlappingBins(location.getContig(),(int)location.getStart(),(int)location.getStop());
    }

    public List<Chunk> getFilePointersBounding(final Bin bin) {
        return reader.getFilePointersBounding(bin);
    }    

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    public int getNumIndexLevels() {
        if(!hasIndex())
            throw new SAMException("Unable to determine number of index levels; BAM file index is not present.");
        return reader.getNumIndexLevels();
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin for which to determine the level.
     * @return the level associated with the given bin number.
     */
    public int getLevelForBin(final Bin bin) {
        if(!hasIndex())
            throw new SAMException("Unable to determine level of bin number; BAM file index is not present.");
        return reader.getLevelForBin(bin);
    }

    public StingSAMIterator seek(Shard shard) {
        if(!(shard instanceof BAMFormatAwareShard))
            throw new StingException("BlockDrivenSAMDataSource cannot operate on shards of type: " + shard.getClass());

        // Since the beginning of time for the GATK, enableVerification has been true only for ReadShards.  I don't
        // know why this is.  Please add a comment here if you do.
        boolean enableVerification = shard instanceof ReadShard;

        CloseableIterator<SAMRecord> iterator = reader.iterator(((BAMFormatAwareShard)shard).getChunks());
        return applyDecoratingIterators(enableVerification,
                StingSAMIteratorAdapter.adapt(reads,iterator),
                reads.getDownsamplingFraction(),
                reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                reads.getSupplementalFilters());
    }

    /**
     * Gets the merged header from the SAM file.
     * @return The merged header.
     */
    public SAMFileHeader getHeader() {
        return reader.getFileHeader();
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
