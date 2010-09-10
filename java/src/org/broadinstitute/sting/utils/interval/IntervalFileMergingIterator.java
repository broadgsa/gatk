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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.interval;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.exceptions.UserError;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.refdata.utils.StringToGenomeLocIteratorAdapter;

import java.util.Iterator;
import java.io.File;
import java.io.FileNotFoundException;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jun 11, 2010
 * Time: 2:56:29 PM
 * To change this template use File | Settings | File Templates.
 */

/** This iterator reads intervals from interval file (can be gatk-style
 * interval list or a bed file) and merges them on the fly. Very much alike
 * IntervalUtils.sortAndMergeIntervals() but the list is read sequentially
 * from a file upon request instead of loading the whole list into memory.
 * Intervals in the underlying file MUST be
 * pre-sorted into the reference order (they can overlap though, as this
 * iterator is a merging one).
 */
public class IntervalFileMergingIterator implements Iterator<GenomeLoc> {
    private PushbackIterator<GenomeLoc> it ;
    private IntervalMergingRule myRule;
    private File myFile;

    public IntervalFileMergingIterator(File f, IntervalMergingRule rule) {
        myFile = f;

        try {
            XReadLines reader = new XReadLines(f);

            if (f.getName().toUpperCase().endsWith(".BED")) {
                it = new PushbackIterator<GenomeLoc>( new StringToGenomeLocIteratorAdapter( reader.iterator(),
                                                              StringToGenomeLocIteratorAdapter.FORMAT.BED ) ) ;
            } else {
                it = new PushbackIterator<GenomeLoc>( new StringToGenomeLocIteratorAdapter( reader.iterator(),
                                                              StringToGenomeLocIteratorAdapter.FORMAT.GATK ) ) ;
            }
        } catch ( FileNotFoundException e ) {
            throw new UserError.CouldNotReadInputFile(f, e);
        }
        myRule = rule;
    }

    public boolean hasNext() {
        return it.hasNext();
    }

    /** Returns next merged interval from the underlying interval file. In other words, keeps reading intervals
     * for as long as they overlap and returns a single merged interval encompassing the set of overlapping
     * intervals read from the file. Non-overlapping intervals are returned as is. This method will throw an
     * exception if it runs into an interval that is out of order.
     * @return
     */
    public GenomeLoc next() {

        GenomeLoc current = it.next();

        while ( it.hasNext() ) {
            GenomeLoc next = it.next();

            if ( next.isBefore(current)) {
                throw new UserError.MalformedFile(myFile, "Interval "+next+" in the interval file is out of order.");
            }

            if (current.overlapsP(next)) {
                current = current.merge(next);
            } else if (current.contiguousP(next) && myRule == IntervalMergingRule.ALL) {
                current = current.merge(next);
            } else {
                it.pushback(next);
                break;
            }
        }

        return current;
    }

    public void remove() {
        throw new UnsupportedOperationException("method 'remove' is not supported by this iterator");
    }

}
