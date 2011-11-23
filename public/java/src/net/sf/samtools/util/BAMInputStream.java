/*
 * Copyright (c) 2011, The Broad Institute
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

package net.sf.samtools.util;

import java.io.IOException;

/**
 * An input stream formulated for use reading BAM files.  Supports
 */
public interface BAMInputStream {
    /**
     * Seek to the given position in the file.  Note that pos is a special virtual file pointer,
     * not an actual byte offset.
     *
     * @param pos virtual file pointer
     */
    public void seek(final long pos) throws IOException;

    /**
     * @return virtual file pointer that can be passed to seek() to return to the current position.  This is
     * not an actual byte offset, so arithmetic on file pointers cannot be done to determine the distance between
     * the two.
     */
    public long getFilePointer();

    /**
     * Determines whether or not the inflater will re-calculated the CRC on the decompressed data
     * and check it against the value stored in the GZIP header.  CRC checking is an expensive
     * operation and should be used accordingly.
     */
    public void setCheckCrcs(final boolean check);

    public int read() throws java.io.IOException;

    public int read(byte[] bytes) throws java.io.IOException;

    public int read(byte[] bytes, int i, int i1) throws java.io.IOException;

    public long skip(long l) throws java.io.IOException;

    public int available() throws java.io.IOException;

    public void close() throws java.io.IOException;

    public void mark(int i);

    public void reset() throws java.io.IOException;

    public boolean markSupported();
}
