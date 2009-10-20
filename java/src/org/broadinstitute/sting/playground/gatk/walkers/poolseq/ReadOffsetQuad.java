package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
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
