package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import com.google.java.contract.*;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 4/10/11
 * Time: 8:49 AM
 *
 * A general interface for ReadCompressors.  Read compressors have the following semantics:
 *
 * The accept a stream of reads, in order, and after each added read returns a compressed stream
 * of reads for emission.  This stream of reads is a "reduced" representation of the total stream
 * of reads.  The actual compression approach is left up to the implementing class.
 */
public interface ConsensusReadCompressor {
    /**
     * Adds the read to the compressor.  The returned iteratable collection of
     * reads represents the incremental compressed output.
     * @param read the next uncompressed read in the input stream to the compressor
     * @return an iterator over the incrementally available compressed reads
     */
    @Requires("read != null")
    @Ensures("result != null")
    Iterable<SAMRecord> addAlignment(SAMRecord read);

    /**
     * Must be called after the last read has been added to finalize the compressor state
     * and return the last compressed reads from the compressor.
     * @return an iterator over the final compressed reads of this compressor
     */
    @Ensures("result != null")
    Iterable<SAMRecord> close();
}
