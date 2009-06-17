package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;

public class RecalData {
    long N;
    long B;
    int pos;
    int qual;
    String readGroup;
    String dinuc;

    public RecalData(int pos, int qual, String readGroup, String dinuc ) {
        this.pos = pos;
        this.qual = qual;
        this.readGroup = readGroup;
        this.dinuc = dinuc;
    }

    public void inc(long incN, long incB) {
        N += incN;
        B += incB;
    }

    public int getDinucIndex() {
        return string2dinucIndex(this.dinuc);
    }

    public void inc(char curBase, char ref) {
        inc(1, nuc2num[curBase] == nuc2num[ref] ? 0 : 1);
        //out.printf("%s %s\n", curBase, ref);
    }

    public static String headerString() {
        return ("pos, rg, dinuc, qual, emp_qual, qual_diff, n, b");
    }

    public double empiricalQualDouble() {
        double empiricalQual = -10 * Math.log10((double)B / N);
        if (empiricalQual > QualityUtils.MAX_QUAL_SCORE) empiricalQual = QualityUtils.MAX_QUAL_SCORE;
        return empiricalQual;
    }

    public byte empiricalQualByte() {
        return QualityUtils.probToQual(1.0 - (double)B / N);
    }

    public String toString() {
        double empiricalQual = empiricalQualDouble();
        return String.format("%3d,%s,%s,%3d,%5.1f,%5.1f,%6d,%6d", pos, readGroup, dinuc, qual, empiricalQual, qual-empiricalQual, N, B);
    }

    public String toCSVString() {
        return String.format("%s,%s,%d,%d,%d,%d", readGroup, dinuc, qual, pos, N, B);
    }

    public static RecalData fromCSVString(String s) {
        String[] vals = s.split(",");
        String rg = vals[0];
        String dinuc = vals[1];
        int qual = Integer.parseInt(vals[2]);
        int pos = Integer.parseInt(vals[3]);
        int N = Integer.parseInt(vals[4]);
        int B = Integer.parseInt(vals[5]);
        RecalData datum = new RecalData(pos, qual, rg, dinuc);
        datum.B = B;
        datum.N = N;

        //if ( datum.N > 0 ) System.out.printf("Parsing line [%s] => [%s]%n", s, datum);        

        return datum;
    }

    public static int bases2dinucIndex(char prevBase, char base, boolean Complement) {
        if (!Complement) {
            return nuc2num[prevBase] * 4 + nuc2num[base];
        }else{
            return (3 - nuc2num[prevBase]) * 4 + (3 - nuc2num[base]);
        }
    }

    public static String dinucIndex2bases(int index) {
        char data[] = {num2nuc[index / 4], num2nuc[index % 4]};
        return new String( data );
    }

    public static int string2dinucIndex(String s) {
        return bases2dinucIndex(s.charAt(0), s.charAt(1), false);
    }

    private static int nuc2num[];
    private static char num2nuc[];
    
    static {
        nuc2num = new int[128];
        nuc2num['A'] = 0;
        nuc2num['C'] = 1;
        nuc2num['G'] = 2;
        nuc2num['T'] = 3;
        nuc2num['a'] = 0;
        nuc2num['c'] = 1;
        nuc2num['g'] = 2;
        nuc2num['t'] = 3;

        num2nuc = new char[4];
        num2nuc[0] = 'A';
        num2nuc[1] = 'C';
        num2nuc[2] = 'G';
        num2nuc[3] = 'T';
    }
}