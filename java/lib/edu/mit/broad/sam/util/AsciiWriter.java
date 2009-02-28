/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.util;

import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;

/**
 * Fast (I hope) Writer that converts char to byte merely by casting, rather than charset conversion.
 */
public class AsciiWriter extends Writer {

    private final OutputStream os;
    // Buffer size has not been tuned.
    private final byte[] buffer = new byte[10000];
    private int numBytes;

    public AsciiWriter(final OutputStream os) {
        this.os = os;
        numBytes = 0;
    }

    public void close() throws IOException {
        flush();
        os.close();
    }

    public void flush() throws IOException {
        os.write(buffer, 0, numBytes);
        numBytes = 0;
        os.flush();
    }

    public void write(final char[] chars, int offset, int length) throws IOException {
        while (length > 0) {
            final int charsToConvert = Math.min(length, buffer.length - numBytes);
            StringUtil.charsToBytes(chars, offset, charsToConvert, buffer, numBytes);
            numBytes += charsToConvert;
            offset += charsToConvert;
            length -= charsToConvert;
            if (numBytes == buffer.length) {
                os.write(buffer, 0, numBytes);
                numBytes = 0;
            }
        }
    }
}
