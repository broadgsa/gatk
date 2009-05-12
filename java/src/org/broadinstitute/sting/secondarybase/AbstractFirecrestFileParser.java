/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
//package edu.mit.broad.picard.illumina;
package org.broadinstitute.sting.secondarybase;

import edu.mit.broad.picard.util.BasicTextFileParser;

import java.io.Closeable;
import java.io.File;
import java.io.FilenameFilter;
import java.util.*;

import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.secondarybase.FirecrestReadData;
import org.broadinstitute.sting.secondarybase.FirecrestFilenameComparator;

/**
 * Abstract base class for implementing parsers for various versions of Firecrest output
 */
public abstract class AbstractFirecrestFileParser implements Iterator<FirecrestReadData>, Iterable<FirecrestReadData>, Closeable {
    protected final int lane;
    protected final File firecrestDirectory;
    private FirecrestReadData next = null;
    private boolean iterating = false;

    /**
     * Examine the bustard directory to see if it is valid, and prepare for parsing
     */
    public AbstractFirecrestFileParser(final File firecrestDirectory, final int lane) {
        this.lane = lane;
        this.firecrestDirectory = firecrestDirectory;
    }

    /**
     * @return true if the given bustard directory contains the appropriate files, or at least enough
     * of them so that it appears to be a Firecrest directory corresponding to the version of the concrete
     * FirecrestFileParser implementation.
     */
    public abstract boolean isValidFirecrestDirectory();

    /**
     * Called before iteration begins.  If this method is called when isValidFirecrestDirectory() had
     * return false, it will generate exceptions that may help the user diagnose the problem.
     */
    protected abstract void prepareToIterate();

    /**
     * @return the next read
     */
    protected abstract FirecrestReadData readNext();


    /**
     * @return an iterator over a set of elements of type FirecrestReadData
     */
    public Iterator<FirecrestReadData> iterator() {
        if (iterating) {
            throw new IllegalStateException("iterator() method can only be called once, before the first call to hasNext()");
        }
        prepareToIterate();
        next = readNext();
        iterating = true;
        return this;
    }

    /**
     * @return  true if the iteration has more elements.  Otherwise returns false.
     */
    public boolean hasNext() {
        if (!iterating) {
            iterator();
        }
        return next != null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return  the next element in the iteration
     * @throws java.util.NoSuchElementException
     */
    public FirecrestReadData next() {

        if (!hasNext()) {
            throw new NoSuchElementException("Iteration has no more elements.");
        }

        final FirecrestReadData result = next;
        next = readNext();
        return result;
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
     * Override, e.g. to close parser
     */
    public void close() {
    }

    public int getLane() { return this.lane; }

    /**
     * Convenience method to create a parser for a list of files of the same format that should
     * be parsed in order defined by FirecrestFilenameComparator
     * @param files to be iterated, in arbitrary order
     * @return parser that iterates through the files in the appropriate order
     */
    protected BasicTextFileParser makeParserForTextFiles(final boolean treatGroupedDelimitersAsOne, File[] files) {
        final SortedSet<File> sortedRead1 = new TreeSet<File>(new FirecrestFilenameComparator());
        sortedRead1.addAll(Arrays.asList(files));
        files = sortedRead1.toArray(files);
        return new BasicTextFileParser(treatGroupedDelimitersAsOne, files);
    }

    protected File[] getFilesMatchingRegexp(final String regexp) {
        return firecrestDirectory.listFiles( new FilenameFilter() {
            public boolean accept(final File dir, final String name) {
                return name.matches(regexp);
            }
        });
    }
}
