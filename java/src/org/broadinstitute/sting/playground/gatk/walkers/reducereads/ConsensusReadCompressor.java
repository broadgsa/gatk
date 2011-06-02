package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

import com.google.java.contract.*;
import net.sf.samtools.SAMRecord;

/*
 * Copyright (c) 2009 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

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
