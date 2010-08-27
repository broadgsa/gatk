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

package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import net.sf.samtools.SAMRecord;
import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 12, 2009
 * Time: 1:21:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReadOffsetQuad extends Quad<List<SAMRecord>,List<Integer>,List<SAMRecord>, List<Integer>> {
    /*
     * ReadOffsetQuad separates the user from specifying the types of objects required to
     * store two sets of read/offset pairs in a quad.
     * 
     * Implements methods that return read/offset pairs as AlignmentContexts
     * and allows ReadOffsetQuad to be constructed from two AlignmentContexts
     *
     */

    // constructor that IntelliJ wants
    public ReadOffsetQuad(List<SAMRecord> a, List<Integer> b, List<SAMRecord> c, List<Integer> d) {
        super(a,b,c,d);
    }

    // another constructor that IntelliJ wants
    public ReadOffsetQuad(Pair<List<SAMRecord>,List<Integer>> a, Pair<List<SAMRecord>,List<Integer>> b) {
        super(a,b);
    }

    public ReadOffsetQuad(AlignmentContext a, AlignmentContext b) {
        first = a.getReads();
        second = a.getOffsets();
        third = b.getReads();
        fourth = b.getOffsets();
    }

    public int numReads() {
        return first.size() + third.size();
    }

    public int numReadsFirst() {
        return first.size();
    }

    public int numReadsSecond() {
        return third.size();
    }

    public List<SAMRecord> getFirstReads() {
        return this.first;
    }

    public List<SAMRecord> getSecondReads() {
        return this.third;
    }

    public List<SAMRecord> getReadsCombined() {
        ArrayList<SAMRecord> combined = new ArrayList<SAMRecord>(first);
        combined.addAll(third);
        return combined;
    }

    public List<Integer> getOffsetsCombined() {
        ArrayList<Integer> combined = new ArrayList<Integer>(second);
        combined.addAll(fourth);
        return combined;
    }

    public List<Integer> getFirstOffsets() {
        return this.second;
    }

    public List<Integer> getSecondOffsets() {
        return this.fourth;
    }

    public Pair<List<SAMRecord>,List<Integer>> getFirstReadOffsetPair() {
        return new Pair<List<SAMRecord>,List<Integer>>(first, second);
    }

    public Pair<List<SAMRecord>,List<Integer>> getSecondReadOffsetPair() {
        return new Pair<List<SAMRecord>,List<Integer>>(third,fourth);
    }

    public AlignmentContext getFirstPairAsAlignmentContext(GenomeLoc loc) {
        return new AlignmentContext(loc, first, second);
    }

    public AlignmentContext getSecondPairAsAlignmentContext(GenomeLoc loc) {
        return new AlignmentContext(loc, third, fourth);
    }

}
