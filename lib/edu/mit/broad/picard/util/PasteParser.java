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

/**
 * Class to merge files horizontally (like the Unix paste command), so that the first line of each file
 * is merged together in one big line, then the second lines, etc.
 *
 * @author Kathleen Tibbetts
 */
public class PasteParser implements Iterable<String[][]>, CloseableIterator<String[][]>{

    private final CloseableIterator<String[]>[] iterators;
    private boolean iterating = false;
    private String[][] next = null;

    /**
     * Constructor
     *
     * @param iterators The iterators containing the files to merge together
     */
    public PasteParser(CloseableIterator<String[]>... iterators) {
        this.iterators = iterators;
    }

    /**
     * Merges the "next" line from each of the underying iterators and returns an array of the results.
     *
     * @return  An array of the lines from each iterator
     * @throws PicardException if the files are not exhausted at the same time
     */
    protected String[][] readNextLine() {
        String result[][] = new String[iterators.length][];
        boolean oneFinished = false;
        boolean oneNotFinished = false;

        for (int i = 0; i < iterators.length; i++) {
            if (!iterators[i].hasNext()) {
                oneFinished = true;
            }
            else {
                result[i] = iterators[i].next();
                oneNotFinished = true;
            }
        }
        if (oneFinished) {
            if (oneNotFinished) {
                throw new PicardException("Mismatched file lengths in PasteParser");
            }
            else {
                return null;
            }
        }
        return result;
    }

    /**
     * Closes the underlying iterators.
     */
    public void close() {
        for (CloseableIterator iterator : iterators) {
            iterator.close();
        }
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
     * Returns an iterator over a set of elements of type BustardReadData.
     *
     * @return an iterator over a set of elements of type BustardReadData
     */
    public Iterator<String[][]> iterator() {
        if (iterating) {
            throw new IllegalStateException("iterator() method can only be called once, before the" +
                    "first call to hasNext()");
        }
        next = readNextLine();
        iterating = true;
        return this;
    }

    /**
     * Returns true if the iteration has more elements.
     *
     * @return  true if the iteration has more elements.  Otherwise returns false.
     */
    public boolean hasNext() {
        if (!iterating) {
            next = readNextLine();
            iterating = true;
        }
        return next != null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return  the next element in the iteration
     * @throws java.util.NoSuchElementException
     */
    public String[][] next() {

        if (!hasNext()) {
            throw new NoSuchElementException("Iteration has no more elements.");
        }

        String[][] result = next;
        next = readNextLine();
        return result;
    }
}
