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

package org.broadinstitute.sting.gatk.contexts;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.*;

/**
 * Useful class for forwarding on locusContext data from this iterator
 * 
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:01:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentContext {
    private GenomeLoc loc = null;
    private List<SAMRecord> reads = null;
    private List<Integer> offsets = null;

    /**
     * Create a new AlignmentContext object
     *
     * @param loc
     * @param reads
     * @param offsets
     */
    public AlignmentContext(GenomeLoc loc, List<SAMRecord> reads, List<Integer> offsets) {
        //assert loc != null;
        //assert loc.getContig() != null;
        //assert reads != null;
        //assert offsets != null;

        this.loc = loc;
        this.reads = reads;
        this.offsets = offsets;
    }

    /**
     * get all of the reads within this context
     * 
     * @return
     */
    public List<SAMRecord> getReads() { return reads; }

    /**
     * Are there any reads associated with this locus?
     *
     * @return
     */
    public boolean hasReads() {
        return reads != null;
    }

    /**
     * How many reads cover this locus?
     * @return
     */
    public int numReads() {
        assert( reads != null );
        return reads.size();
    }

    /**
     * get a list of the equivalent positions within in the reads at Pos
     *
     * @return
     */
    public List<Integer> getOffsets() {
        return offsets;
    }

    public String getContig() { return getLocation().getContig(); }
    public long getPosition() { return getLocation().getStart(); }
    public GenomeLoc getLocation() { return loc; }
    public void setLocation(GenomeLoc loc) {
        this.loc = loc.clone();
    }

    public void downsampleToCoverage(int coverage) {
        if ( numReads() <= coverage )
            return;

        // randomly choose numbers corresponding to positions in the reads list
        Random generator = new Random();
        TreeSet positions = new TreeSet();
        int i = 0;
        while ( i < coverage ) {
            if (positions.add(new Integer(generator.nextInt(reads.size()))))
                i++;
        }

        ArrayList<SAMRecord> downsampledReads = new ArrayList<SAMRecord>();
        ArrayList<Integer> downsampledOffsets = new ArrayList<Integer>();
        Iterator positionIter = positions.iterator();
        Iterator<SAMRecord> readsIter = reads.iterator();
        Iterator<Integer> offsetsIter = offsets.iterator();
        int currentRead = 0;
        while ( positionIter.hasNext() ) {
            int nextReadToKeep = (Integer)positionIter.next();

            // fast-forward to the right read
            while ( currentRead < nextReadToKeep ) {
                readsIter.next();
                offsetsIter.next();
                currentRead++;
            }

            downsampledReads.add(readsIter.next());
            downsampledOffsets.add(offsetsIter.next());
            currentRead++;
        }

        reads = downsampledReads;
        offsets = downsampledOffsets;
    }
}
