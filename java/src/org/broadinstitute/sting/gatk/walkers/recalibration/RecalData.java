package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.List;
import java.util.Collections;

public class RecalData implements Comparable<RecalData> {
    //
    // context of this data -- read group, position in read, reported quality, and dinuc
    //
    String readGroup;       // the read group where the datum was observed
    int pos;                // the position of this recalibration datum
    int qual;               // the reported quality scoree of this datum
    String dinuc;           // the dinucleotide context of this datum

    //
    // empirical quality score data
    //
    long N;                 // total number of observed bases
    long B;                 // number of observed mismatches to the reference base

    
    public RecalData(int pos, int qual, String readGroup, String dinuc ) {
        this.pos = pos;
        this.qual = qual;
        this.readGroup = readGroup;
        this.dinuc = dinuc;
    }

    public RecalData( RecalData copy ) {
        this.pos = copy.pos;
        this.qual = copy.qual;
        this.readGroup = copy.readGroup;
        this.dinuc = copy.dinuc;
        this.N = copy.N;
        this.B = copy.B;
    }

    public int compareTo(RecalData to) {
        int cmpReadGroup = readGroup.compareTo(to.readGroup);
        if (cmpReadGroup != 0) return cmpReadGroup;

        int cmpPos = new Integer(pos).compareTo(to.pos);
        if (cmpPos != 0) return cmpPos;

        int cmpQual = new Integer(qual).compareTo(to.qual);
        if (cmpQual != 0) return cmpQual;

        return dinuc.compareTo(to.dinuc);
    }

    /**
     * Returns the expected number of mismatching bases for this datum, using the number of bases N and
     * the reported qual score.  Effectively qual implied error rate * number of bases.
     * 
     * @return
     */
    public double getExpectedErrors() {
        return QualityUtils.qualToErrorProb((byte)qual) * N;
    }

    //
    // Increment routines
    //
    public RecalData inc(long incN, long incB) {
        N += incN;
        B += incB;
        return this;
    }

    /**
     * Add N and B from other to self
     *
     * @param other
     * @return
     */
    public RecalData inc(final RecalData other) {
        return inc(other.N, other.B);
    }

    /**
     * Add all of the Ns and Bs for data to self
     *
     * @param data
     * @return
     */
    public RecalData inc(List<RecalData> data) {
        for ( RecalData other : data ) {
            this.inc(other);
        }
        return this;
    }

    /**
     * Count an empirical observation of our current base against the reference base.
     * Increments N and, if curBase != ref, also increments B
     *
     * @param curBase
     * @param ref
     * @return
     */
    public RecalData inc(char curBase, char ref) {
        return inc(1, BaseUtils.simpleBaseToBaseIndex(curBase) == BaseUtils.simpleBaseToBaseIndex(ref) ? 0 : 1);
        //out.printf("%s %s\n", curBase, ref);
    }


    /**
     * Get the integer dinuc index associated with this datum.  Can return -1 if dinuc
     * contains unexpected bases.  Returns 0 if dinuc contains the * operator (not tracking dinucs)
     * @return
     */
    public int getDinucIndex() {
        return dinucIndex(this.dinuc);
    }

    /**
     * Calculates the empirical quality score for this datum from B and N.  Bounds the result
     * by QualityUtils.MAX_REASONABLE_Q_SCORE.
     *
     * @return
     */
    public double empiricalQualDouble(final boolean useRawQempirical) {
        double doubleB = useRawQempirical ? B : B + 1;
        double doubleN = useRawQempirical ? N : N + 1;
        double empiricalQual = -10 * Math.log10(doubleB / doubleN);
        if (empiricalQual > QualityUtils.MAX_REASONABLE_Q_SCORE) empiricalQual = QualityUtils.MAX_REASONABLE_Q_SCORE;
        return empiricalQual;
    }
    public double empiricalQualDouble() { return empiricalQualDouble(true); }


    /**
     * As as empiricalQualDouble, but rounded and encoded as a byte for placement in a read quality array of bytes
     * @return
     */
    public byte empiricalQualByte(final boolean useRawQempirical) {
        double doubleB = useRawQempirical ? B : B + 1;
        double doubleN = useRawQempirical ? N : N + 1;
        return QualityUtils.probToQual(1.0 - doubleB / doubleN);
    }
    public byte empiricalQualByte() { return empiricalQualByte(true); }


    public static String headerString() {
        return ("pos, rg, dinuc, qual, emp_qual, qual_diff, n, b");
    }

    public String toString() {
        return String.format("[rg=%s, pos=%3d, Qrep=%3d, dinuc=%s, Qemp=%2.2f / %2.2f, B=%d, N=%d]", readGroup, pos, qual, dinuc, empiricalQualDouble(), empiricalQualDouble(false), B, N);
        //return String.format("%3d,%s,%s,%3d,%5.1f,%5.1f,%6d,%6d", pos, readGroup, dinuc, qual, empiricalQual, qual-empiricalQual, N, B);
    }

    //
    // To and from CSV strings
    //
    public String toCSVString(boolean collapsedPos) {
        // rg,pos,Qrep,dn,nBases,nMismatches,Qemp
        return String.format("%s,%s,%d,%s,%d,%d,%d", readGroup, collapsedPos ? "*" : pos, qual, dinuc, N, B, empiricalQualByte());
    }

    /**
     * Parses the string s and returns the encoded RecalData in it.  Assumes data has been emitted by
     * toCSVString() routine.  Order is readGroup, dinuc, qual, pos, N, B, Qemp as byte
     * @param s
     * @return
     */
    public static RecalData fromCSVString(String s) {
        String[] vals = s.split(",");
        String rg = vals[0];
        int pos = vals[1].equals("*") ? 0 : Integer.parseInt(vals[1]);
        int qual = Integer.parseInt(vals[2]);
        String dinuc = vals[3];
        int N = Integer.parseInt(vals[4]);
        int B = Integer.parseInt(vals[5]);
        RecalData datum = new RecalData(pos, qual, rg, dinuc);
        datum.B = B;
        datum.N = N;
        return datum;
    }


    public static List<RecalData> sort(List<RecalData> l) {
        Collections.sort(l);
        return l;
    }

    public static int bases2dinucIndex(char prevBase, char base) {
        int pbI = BaseUtils.simpleBaseToBaseIndex(prevBase);
        int bI = BaseUtils.simpleBaseToBaseIndex(base);
        return (pbI == -1 || bI == -1) ? -1 : pbI * 4 + bI;
    }

    public final static int NDINUCS = 16;
    public static String dinucIndex2bases(int index) {
        char data[] = {BaseUtils.baseIndexToSimpleBase(index / 4), BaseUtils.baseIndexToSimpleBase(index % 4)};
        return new String( data );
    }

    public static int dinucIndex(String s) {
        return bases2dinucIndex(s.charAt(0), s.charAt(1));
    }

    public static int dinucIndex(byte prevBase, byte base) {
        return bases2dinucIndex((char)prevBase, (char)base);
    }

    public static double combinedQreported(List<RecalData> l) {
        double sumExpectedErrors = 0;
        long nBases = 0;
        for ( RecalData datum : l ) {
            sumExpectedErrors += datum.getExpectedErrors();
            nBases += datum.N;
        }

        double q = QualityUtils.phredScaleErrorRate(sumExpectedErrors / nBases);
        System.out.printf("expected errors=%f, nBases = %d, rate=%f, qual=%f%n",
                sumExpectedErrors, nBases, 1 - sumExpectedErrors / nBases, q);
        return q;
    }
}