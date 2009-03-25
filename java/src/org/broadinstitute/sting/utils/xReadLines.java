package org.broadinstitute.sting.utils;

import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Mar 25, 2009
 * Time: 10:46:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class xReadLines implements Iterator<String>, Iterable<String> {
    BufferedReader in;          // The stream we're reading from
    String nextline = null;     // Return value of next call to next()

    public xReadLines(final File filename) throws FileNotFoundException {
        // Open the file and read and remember the first line.
        // We peek ahead like this for the benefit of hasNext().
        this(new FileReader(filename));
    }

    public xReadLines(final FileReader fileReader) throws FileNotFoundException {
        // Open the file and read and remember the first line.
        // We peek ahead like this for the benefit of hasNext().
        this(new BufferedReader(fileReader));
    }

    public xReadLines(final InputStream inputStream) throws FileNotFoundException {
        // Open the file and read and remember the first line.
        // We peek ahead like this for the benefit of hasNext().
        this(new BufferedReader(new InputStreamReader(inputStream)));
    }

    public xReadLines(final BufferedReader in) throws FileNotFoundException {
        // Open the file and read and remember the first line.
        // We peek ahead like this for the benefit of hasNext().
        try {
            this.in = in;
            nextline = readNextLine();
        } catch(IOException e) {
            throw new IllegalArgumentException(e);
        }
    }

    public List<String> readLines() {
        List<String> lines = new LinkedList<String>();
        for ( String line : this ) {
            lines.add(line);
        }
        return lines;
    }

    public Iterator<String> iterator() {
        return this;
    }

    // If the next line is non-null, then we have a next line
    public boolean hasNext() {
        return nextline != null;
    }

    private String readNextLine() throws IOException {
        String nextline = in.readLine();   // Read another line
        if (nextline == null)
            return null;
        else
            nextline = nextline.trim();
        return nextline;
    }

    // Return the next line, but first read the line that follows it.
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
