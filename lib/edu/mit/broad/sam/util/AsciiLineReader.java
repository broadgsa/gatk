/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This is copyright (2007-2008) by the Broad Institute/Massachusetts Institute
 * of Technology.  It is licensed to You under the Gnu Public License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 *  the License.  You may obtain a copy of the License at
 *
 *    http://www.opensource.org/licenses/gpl-2.0.php
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever. Neither the Broad Institute nor MIT can be responsible for its
 * use, misuse, or functionality.
 */
package edu.mit.broad.sam.util;

import java.io.IOException;
import java.io.InputStream;

/**
 * Fast replacement for BufferedReader that assumes that bytes can be converted to chars simply by casting.
 * @author jrobinso
 */
public class AsciiLineReader implements LineReader {
    private static final byte LINEFEED = (byte)('\n' & 0xff);
    private static final byte CARRIAGE_RETURN = (byte)('\r' & 0xff);

    private final InputStream is;
    private byte[] buffer;
    private int nextChar;
    private int nChars;
    // Allocate this only once, despite the fact that it is essentially a local variable of readLine()
    private byte[] lineBuffer = new byte[1000];

    private int lineNumber = 0;

    public AsciiLineReader(final InputStream is) {
        this(is, 512000);
    }

    public AsciiLineReader(final InputStream is, final int bufferSize) {
        this.is = is;
        buffer = new byte[bufferSize];
        nextChar = nChars = 0;
    }

    public String readLine() {
        return readLine(false);
    }

    /**
     * Read a line of text.  A line is considered to be terminated by any one
     * of a line feed ('\n'), a carriage return ('\r'), or a carriage return
     * followed immediately by a linefeed.
     *
     * @param      includeTerminators  If true, the line-termination characters
     *             are included in the returned string. 
     *
     * @return     A String containing the contents of the line or null if the 
     *             end of the stream has been reached
     */
    public String readLine(final boolean includeTerminators){
        int linePosition = 0;

        while (true)
        {
            if (nChars == -1)
            {
                return null;
            }

            // Refill buffer if neccessary
            if (nextChar == nChars)
            {
                fill();
                if (nextChar == nChars || nChars == -1)
                {
                    // eof reached.  Return the last line, or null if this is a new line
                    if (linePosition > 0)
                    {
                        ++lineNumber;
                        return StringUtil.bytesToString(lineBuffer, 0, linePosition);
                    } else
                    {
                        return null;
                    }
                }
            }


            final byte b = buffer[nextChar++];
            if (b == LINEFEED || b == CARRIAGE_RETURN)
            {

                if (includeTerminators)
                {
                    lineBuffer[linePosition++] = b;
                    if (b == CARRIAGE_RETURN && peek() == LINEFEED)
                    {
                        lineBuffer[linePosition++] = b;
                        nextChar++; // <= to account for the '\n' we just ate
                    }
                }
                else {
                    if (b == CARRIAGE_RETURN && peek() == LINEFEED)
                    {
                        nextChar++; // <= skip the trailing \n in case of \r\n termination
                    }
                    
                }
                ++lineNumber;
                return StringUtil.bytesToString(lineBuffer, 0, linePosition);
            } else
            {
                // Expand line buffer size if neccessary.  Reservce at least 2 characters
                // for potential line-terminators in return string

                if (linePosition > (lineBuffer.length - 3))
                {
                    final byte[] temp = new byte[lineBuffer.length + 100];
                    System.arraycopy(lineBuffer, 0, temp, 0, lineBuffer.length);
                    lineBuffer = temp;
                }

                lineBuffer[linePosition++] = b;
            }
        }
    }

    public int getLineNumber() {
        return lineNumber;
    }

    /**
     * Peek ahead one character, filling from the underlying stream if neccessary.
     * 
     * @return
     * @throws java.io.IOException
     */
    private byte peek(){
        // Refill buffer if neccessary
        if (nextChar == nChars)
        {
            fill();
            if (nextChar == nChars)
            {
                // eof reached.  
                return 0;
            }
        }
        return buffer[nextChar];

    }

    private void fill() {
        try {
            nChars = is.read(buffer);
            nextChar = 0;
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }
    
    public void close()  {
        try {
            is.close();
        } catch (IOException e) {
            // Ignore exception
        }
    }
}

