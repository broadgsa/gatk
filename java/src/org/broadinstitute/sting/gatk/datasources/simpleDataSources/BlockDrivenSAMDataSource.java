package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.BlockDelimitedReadShard;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.utils.StingException;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader2;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.util.Collection;
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

        if(reads.getReadsFiles().size() > 1)
            throw new StingException("Experimental sharding strategy cannot handle multiple BAM files at this point.");

        File readsFile = reads.getReadsFiles().get(0);
        reader = new SAMFileReader2(readsFile);
    }

    public boolean hasIndex() {
        return reader.hasIndex();
    }

    public StingSAMIterator seek(Shard shard) {
        if(!(shard instanceof BlockDelimitedReadShard))
            throw new StingException("Currently unable to operate on types other than block delimited read shards.");
        CloseableIterator<SAMRecord> iterator = reader.iterator(((BlockDelimitedReadShard)shard).getChunks());
        return applyDecoratingIterators(true,
                                        StingSAMIteratorAdapter.adapt(reads, iterator),
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
