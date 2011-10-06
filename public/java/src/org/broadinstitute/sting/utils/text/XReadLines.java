/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.text;

import java.io.*;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * Support for Python-like xreadlines() function as a class.  This is an iterator and iterable over
 * Strings, each corresponding a line in the file (minus newline).  Enables the very simple accessing
 * of lines in a file as:
 *
 * xReadLines reader = new xReadLines(new File(file_name));
 * List<String> lines = reader.readLines();
 * reader.close();
 *
 * or
 *
 * for ( String line : new xReadLines(new File(file_name)) {
 *   doSomeWork(line);
 * }
 *
 * For the love of god, please use this system for reading lines in a file.
 */
public class XReadLines implements Iterator<String>, Iterable<String> {
    private BufferedReader in;          // The stream we're reading from
    private String nextline = null;     // Return value of next call to next()
    private boolean trimWhitespace = true;

    /**
     * Creates a new xReadLines object to read lines from filename
     *
     * @param filename
     * @throws FileNotFoundException
     */
    public XReadLines(final File filename, final boolean trimWhitespace) throws FileNotFoundException {
        this(new FileReader(filename), trimWhitespace);
    }

    public XReadLines(final File filename) throws FileNotFoundException {
        this(filename, true);
    }

    /**
     * Creates a new xReadLines object to read lines from fileReader
     *
     * @param fileReader
     * @throws FileNotFoundException
     */
    public XReadLines(final FileReader fileReader, final boolean trimWhitespace) throws FileNotFoundException {
        this(new BufferedReader(fileReader), trimWhitespace);
    }

    public XReadLines(final FileReader fileReader) throws FileNotFoundException {
        this(fileReader, true);
    }

    /**
     * Creates a new xReadLines object to read lines from an input stream
     *
     * @param inputStream
     */
    public XReadLines(final InputStream inputStream, final boolean trimWhitespace) {
        this(new BufferedReader(new InputStreamReader(inputStream)), trimWhitespace);
    }

    public XReadLines(final InputStream inputStream) throws FileNotFoundException {
        this(inputStream, true);
    }


    /**
     * Creates a new xReadLines object to read lines from an bufferedReader
     *
     * @param reader
     */
    public XReadLines(final Reader reader, final boolean trimWhitespace) {
        try {
            this.in = new BufferedReader(reader);
            nextline = readNextLine();
            this.trimWhitespace = trimWhitespace;
        } catch(IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    public XReadLines(final Reader reader) {
        this(reader, true);
    }

    /**
     * Reads all of the lines in the file, and returns them as a list of strings
     *
     * @return
     */
    public List<String> readLines() {
        List<String> lines = new LinkedList<String>();
        for ( String line : this ) {
            lines.add(line);
        }
        return lines;
    }

    /**
     * I'm an iterator too...
     * @return
     */
    public Iterator<String> iterator() {
        return this;
    }

    public boolean hasNext() {
        return nextline != null;
    }

    /**
     * Actually reads the next line from the stream, not accessible publically
     * @return
     */
    private String readNextLine() throws IOException {
        String nextline = in.readLine();   // Read another line
        if (nextline != null && trimWhitespace )
            nextline = nextline.trim();
        return nextline;
    }

    /**
     * Returns the next line (minus whitespace) 
     * @return
     */
    public String next() {
        try {
            String result = nextline;
            nextline = readNextLine();

            // If we haven't reached EOF yet
            if (nextline == null) {
                in.close();             // And close on EOF
            }

            // Return the line we read last time through.
            return result;
        } catch(IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    // The file is read-only; we don't allow lines to be removed.
    public void remove() {
        throw new UnsupportedOperationException();
    }

    public void close() throws IOException {
        this.in.close();
    }
}
