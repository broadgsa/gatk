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

package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Writes an occurrence array to the output file.
 *
 * @author mhanna
 * @version 0.1
 */
public class OccurrenceOutputStream {
    /**
     * How many bytes does it take to hold an integer in Java?
     */
    private static final int INT_SIZE_IN_BYTES = 4; 

    /**
     * Ultimate target for the occurrence array.
     */
    private final OutputStream targetOutputStream;

    /**
     * Create a new OccurrenceArrayOutputStream, writing to the given target file.
     * @param outputFile target output file.
     * @throws IOException if an I/O error occurs.
     */
    public OccurrenceOutputStream( File outputFile ) throws IOException {
        this(new FileOutputStream(outputFile));
    }

    /**
     * Write occurrence array to the given OutputStream.
     * @param outputStream Output stream to which to write packed bases.
     * @throws IOException if an I/O error occurs.
     */
    public OccurrenceOutputStream( OutputStream outputStream ) throws IOException {
        this.targetOutputStream = outputStream;
    }

    /**
     * Write the cumulative occurrences to the output stream.
     * @param occurrences occurrences to write.  occurrences.length must match alphabet size.
     * @throws IOException if an I/O error occurs.
     */
    public void write( int[] occurrences ) throws IOException {
        if( occurrences.length > BytePackedOutputStream.ALPHABET_SIZE )
            throw new StingException("Wrong number of occurrence data points; expected " + BytePackedOutputStream.ALPHABET_SIZE);
        ByteBuffer buffer = ByteBuffer.allocate(INT_SIZE_IN_BYTES*occurrences.length).order(ByteOrder.LITTLE_ENDIAN);
        for(int occurrence: occurrences)
            buffer.putInt(occurrence);
        targetOutputStream.write(buffer.array());        
    }

    /**
     * Flush the contents of the OutputStream to disk.
     * @throws IOException if an I/O error occurs.
     */
    public void flush() throws IOException {
        targetOutputStream.flush();
    }    

    /**
     * Closes the given output stream.
     * @throws IOException if an I/O error occurs.
     */
    public void close() throws IOException {
        targetOutputStream.close();
    }

}
