/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.util;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.sam.util.CloseableIterator;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.io.Closeable;

/**
 * Class for parsing text files where each line consists of fields separated by whitespace.
 * Code is abstracted into this class so that we can optimize its performance over time.
 *
 * This class assumes that every line will have the same number of whitespace-separated "words"
 * and that lines that start with "#" are comments and should be ignored.
 *
 * Classes that extend this parser can do so simply by implementing their own constructors and the
 * readNextLine(), close(), and getFileName() methods.
 *
 * @author Kathleen Tibbetts
 */
public abstract class AbstractTextFileParser implements Iterable<String[]>, CloseableIterator<String[]> {

    private boolean treatGroupedDelimitersAsOne = true; // Whether multiple delimiters in succession should be treated as one
    private byte nextLine[] = null;
    private int wordCount = 0;      /* The number of delimiter-separated "words" per line of the file.
                                       We can save a little caclulation, or handle files with varying numbers of
                                       words per line, by specifying this if known in advance */
    private boolean iterating = false;

    /**
     * Closes this stream and releases any system resources associated with it.
     */
    public abstract void close();

    /**
     * @return the next line of text from the underlying stream(s) or null if there is no next line
     */
    protected abstract byte[] readNextLine();

    /**
     * @return  the name(s) of the file(s) being parsed, or null if no name is available
     */
    protected abstract String getFileName();

    /**
     * @return an iterator over a set of elements of type String[]
     */
    public Iterator<String[]> iterator() {
        if (iterating) {
            throw new IllegalStateException("iterator() method can only be called once, before the" +
                    "first call to hasNext()");
        }
        nextLine = readNextLine();
        iterating = true;
        return this;
    }

    /**
     * Returns true if the iteration has more elements.
     *
     * @return  true if the iteration has more elements.  Otherwise returns false.
     */
    public boolean hasNext() {
        // If this is the start of iteration, queue up the first item
        if(!iterating) {
            nextLine = readNextLine();
            iterating = true;
        }
        return nextLine != null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return  the next tlement in the iteration
     * @throws java.util.NoSuchElementException
     */
    public String[] next() {

        if (!hasNext()) {
            throw new NoSuchElementException("Iteration from text file(s) " +
                    getFileName() + " has no more elements.");
        }

        String[] result = parseLine(nextLine);
        do {
            nextLine = readNextLine();
        }
        while (nextLine != null && isComment(nextLine));
        return result;
    }

    /**
     * This method represents the most efficient way (so far) to parse a line of whitespace-delimited text
     *
     * @param line the line to parse
     * @return  an array of all the "words"
     */
    private String[] parseLine(byte line[]) {

        if (getWordCount() == 0) {
            calculateWordCount(line);
        }
        String parts[] = new String[getWordCount()];
        boolean delimiter = true;
        int index=0;
        int start = 0;

        try
        {
            for (int i = 0; i < line.length; i++) {
                if (isDelimiter(line[i])) {
                    if (!delimiter) {
                        parts[index++] = new String(line,start,i-start);
                    }
                    else if(!isTreatGroupedDelimitersAsOne()) {
                        parts[index++] = null;
                    }
                    delimiter=true;
                }
                else {
                    if (delimiter)  start = i;
                    delimiter = false;
                }
            }
            if (!delimiter) {
                 parts[index] = new String(line,start,line.length-start);
            }
        }
        catch (ArrayIndexOutOfBoundsException e) {
            throw new PicardException("Unexpected number of elements found when parsing file " +
                    this.getFileName() + ": " + index + ".  Expected a maximum of " +
                    this.getWordCount() + " elements per line.");
        }
        return parts;
    }

    /**
     * Calculates the number of delimiter-separated "words" in a line and sets the value of <code>wordCount</code>
     *
     * @param line  representative line from the file
     */
    protected void calculateWordCount(byte line[]) {
        int words = 0;
        boolean delimiter = true;
        for (byte b : line) {
            if (isDelimiter(b)) {
                if (delimiter && !isTreatGroupedDelimitersAsOne()) words++;
                delimiter = true;
            } else {
                if (delimiter) words++;
                delimiter = false;
            }
        }
        setWordCount(words);
    }

    /**
     * Required method for Iterator API.
     *
     * @throws UnsupportedOperationException
     */
    public void remove() {
        throw new UnsupportedOperationException("Remove() not supported.");
    }

    /**
     * Determines whether a given line is a comment
     *
     * @param line  the line to evaluate
     * @return  true if the line is a comment (and should be ignored) otherwise false
     */
    protected boolean isComment(byte line[]) {
        return line[0] == '#';
    }

    /**
     * Determines whether a given character is a delimiter
     *
     * @param b the character to evaluate
     * @return  true if <code>b</code> is a delimiter; otherwise false
     */
    protected boolean isDelimiter(byte b) { 
        return b == ' ' || b == '\t';
    }

    protected int getWordCount() { return wordCount; }
    protected void setWordCount(int wordCount) { this.wordCount = wordCount; }
    protected boolean isTreatGroupedDelimitersAsOne() { return treatGroupedDelimitersAsOne; }
    protected void setTreatGroupedDelimitersAsOne(boolean treatGroupedDelimitersAsOne) {
        this.treatGroupedDelimitersAsOne = treatGroupedDelimitersAsOne;
    }
}
