package org.broadinstitute.sting.playground.indels;

import org.broadinstitute.sting.utils.PrimitivePair;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.Cigar;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 23, 2009
 * Time: 1:54:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class SWPairwiseAlignment {
    private String s1;
    private String s2;
    private int i1;
    private int i2;
    private int alignment_offset; // offset of s2 w/respect to s1
    private Cigar alignmentCigar;

    private int best_mm; // mismatch count
    private static final int IMPOSSIBLE = 1000000000;
    private static final int MSTATE = 0;
    private static final int ISTATE = 1;
    private static final int DSTATE = 2;

    public SWPairwiseAlignment(String seq1, String seq2, int id1, int id2 ) {
        s1 = seq1;
        s2 = seq2;
        i1 = id1;
        i2 = id2;
        best_mm = IMPOSSIBLE;
        //next_mm = IMPOSSIBLE;
        align2(s1,s2);
    }

    /** Initializes the alignment with pair of sequences (that will be immediately aligned) and
     * sets their external ids to -1. Such un-annotated pairwise alignment can not be added to MultipleAlignment.
     *
     */
    public SWPairwiseAlignment(String seq1, String seq2) {
        this(seq1,seq2,-1,-1);
    }

    public Cigar getCigar() { return alignmentCigar ; }

    public int getAlignmentStart2wrt1() { return alignment_offset; }

    public void align(String a, String b) {
        int n = a.length();
        int m = b.length();
        int [][] sw = new int[n+1][m+1];

        // build smith-waterman matrix:
        for ( int i = 1 ; i < n+1 ; i++ ) {
            char a_base = Character.toUpperCase(a.charAt(i-1)); // letter in a at the current pos
            for ( int j = 1 ; j < m+1 ; j++ ) {
                char b_base = Character.toUpperCase(b.charAt(j-1)); // letter in b at the current pos
                int step_diag = sw[i-1][j-1] + w(a_base,b_base);
                int step_down = sw[i-1][j]+w(a_base,'-');
                int step_right = sw[i][j-1]+w('-',b_base);

                sw[i][j] = Math.max(0, Math.max(step_diag,Math.max(step_down,step_right)));
            }
        }

        print(sw,a,b);

        PrimitivePair.Int p = new PrimitivePair.Int();
        int maxscore = 0;
        // look for largest score. we use >= combined with the traversal direction
        // to ensure that if two scores are equal, the one closer to diagonal gets picked
        for ( int i = 1 ; i < n+1 ; i++ ) {
            if ( sw[i][m] >= maxscore ) { p.first = i; p.second = m ; maxscore = sw[i][m]; }
        }
        for ( int j = 1 ; j < m+1 ; j++ ) {
            if ( sw[n][j] > maxscore ||
                  sw[n][j] == maxscore && Math.abs(n-j) < Math.abs(p.first-p.second)) {
                p.first = n;
                p.second = j ;
                maxscore = sw[n][j];
            }
        }

        System.out.println("\ni="+p.first+"; j="+p.second);

        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order


        // we will be placing all insertions and deletions into sequence b, so the state are named w/regard
        // to that sequence

        int state = MSTATE;
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

        int [] scores = new int[3];

        List<CigarElement> lce = new ArrayList<CigarElement>(5);

        do {
            scores[ISTATE] = sw[p.first][p.second-1]; // moving left: same base on a, prev base on b = insertion on b
            scores[DSTATE] = sw[p.first-1][p.second]; // moving up: same base on b, prev base on a = deletion on b
            scores[MSTATE] = sw[p.first-1][p.second-1]; // moving diagonal : mathc/mismatch

            int new_state = findMaxInd(scores,MSTATE);

            // move to next best location in the sw matrix:
            switch( new_state ) {
                case MSTATE: p.first--; p.second--; break;
                case ISTATE: p.second--; break;
                case DSTATE: p.first--; break;
            }

            // now let's see if the state actually changed:
            if ( new_state == state ) segment_length++;
            else {
                // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
                CigarOperator o=null;
                switch(state) {
                    case MSTATE: o = CigarOperator.M; break;
                    case ISTATE: o = CigarOperator.I; break;
                    case DSTATE: o = CigarOperator.D; break;
                }
                CigarElement e = new CigarElement(segment_length,o);
                lce.add(e);
                segment_length = 1;
                state = new_state;
            }
        } while ( scores[state] != 0 );

        // post-process the last segment we are still keeping
        CigarOperator o=null;

        switch(state) {
            case MSTATE: o = CigarOperator.M; break;
            case ISTATE: o = CigarOperator.I; break;
            case DSTATE: o = CigarOperator.D; break;
        }

        CigarElement e = new CigarElement(segment_length,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
        alignment_offset = p.first - p.second;
    }

    /** Allows for separate gap opening end extension penalties, no direct backtracking.
     *
     * @param a
     * @param b
     */
    public void align2(String a, String b) {
        int n = a.length();
        int m = b.length();
        double [][] sw = new double[n+1][m+1];

        // build smith-waterman matrix:
        for ( int i = 1 ; i < n+1 ; i++ ) {
            char a_base = Character.toUpperCase(a.charAt(i-1)); // letter in a at the current pos
            for ( int j = 1 ; j < m+1 ; j++ ) {
                char b_base = Character.toUpperCase(b.charAt(j-1)); // letter in b at the current pos
                double step_diag = sw[i-1][j-1] + wd(a_base,b_base);
                double step_down = 0.0 ;
                for ( int k = 1 ; k < i ; k++ ) step_down = Math.max(step_down,sw[i-k][j]+wk(a_base,'-',k));
                    
                double step_right = 0;
                for ( int k = 1 ; k < j ; k++ ) step_right = Math.max(step_right,sw[i][j-k]+wk('-',b_base,k));

                sw[i][j] = Math.max(0, Math.max(step_diag,Math.max(step_down,step_right)));
            }
        }

        print(sw,a,b);

        PrimitivePair.Int p = new PrimitivePair.Int();
        double maxscore = 0.0;
        // look for largest score. we use >= combined with the traversal direction
        // to ensure that if two scores are equal, the one closer to diagonal gets picked
        for ( int i = 1 ; i < n+1 ; i++ ) {
            if ( sw[i][m] >= maxscore ) { p.first = i; p.second = m ; maxscore = sw[i][m]; }
        }
        for ( int j = 1 ; j < m+1 ; j++ ) {
            if ( sw[n][j] > maxscore ||
                  sw[n][j] == maxscore && Math.abs(n-j) < Math.abs(p.first-p.second)) {
                p.first = n;
                p.second = j ;
                maxscore = sw[n][j];
            }
        }

        System.out.println("\ni="+p.first+"; j="+p.second);

        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order


        // we will be placing all insertions and deletions into sequence b, so the state are named w/regard
        // to that sequence

        int state = MSTATE;
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

        double [] scores = new double[3];

        List<CigarElement> lce = new ArrayList<CigarElement>(5);

        do {
            // moving left: same base on a, prev base on b = insertion on b:
            scores[ISTATE] = sw[p.first][p.second-1] ;
            scores[DSTATE] = sw[p.first - 1][p.second];
            scores[MSTATE] = sw[p.first-1][p.second-1]; // moving diagonal : match/mismatch

            //            System.out.println("i = " + p.first + " ; j = " + p.second);
            //            System.out.println("s(M)="+scores[MSTATE]+"; s(D)="+scores[DSTATE]+"; s(I)=" + scores[ISTATE]);
            int new_state = findMaxInd(scores,MSTATE);

            // move to next best location in the sw matrix:
            switch( new_state ) {
                case MSTATE: p.first--; p.second--; break;
                case ISTATE: p.second--; break;
                case DSTATE: p.first--; break;
            }

            // now let's see if the state actually changed:
            if ( new_state == state ) segment_length++;
            else {
                // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
                CigarOperator o=null;
                switch(state) {
                    case MSTATE: o = CigarOperator.M; break;
                    case ISTATE: o = CigarOperator.I; break;
                    case DSTATE: o = CigarOperator.D; break;
                }
                CigarElement e = new CigarElement(segment_length,o);
                lce.add(e);
                segment_length = 1;
                state = new_state;
            }
        } while ( scores[state] != 0 );

        // post-process the last segment we are still keeping
        CigarOperator o=null;
        switch(state) {
            case MSTATE: o = CigarOperator.M; break;
            case ISTATE: o = CigarOperator.I; break;
            case DSTATE: o = CigarOperator.D; break;
        }
        CigarElement e = new CigarElement(segment_length,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
        alignment_offset = p.first - p.second ;
    }


    /** Allows for separate gap opening and extension penalties, with backtracking.
     *
     * @param a
     * @param b
     */
public void align3(String a, String b) {
        int n = a.length();
        int m = b.length();
        double [][] sw = new double[n+1][m+1];

        int [][] btrack = new int[n+1][m+1];

        // build smith-waterman matrix and keep backtrack info:
        for ( int i = 1 ; i < n+1 ; i++ ) {
            char a_base = Character.toUpperCase(a.charAt(i-1)); // letter in a at the current pos
            for ( int j = 1 ; j < m+1 ; j++ ) {
                char b_base = Character.toUpperCase(b.charAt(j-1)); // letter in b at the current pos
                double step_diag = sw[i-1][j-1] + wd(a_base,b_base);
                double step_down = 0.0 ;
                int kd = 0;
                for ( int k = 1 ; k < i ; k++ ) {
                    if ( step_down < sw[i-k][j]+wk(a_base,'-',k) ) {
                        step_down=sw[i-k][j]+wk(a_base,'-',k);
                        kd = k;
                    }
                }

                double step_right = 0;
                int ki = 0;
                for ( int k = 1 ; k < j ; k++ ) {
                    if ( step_right < sw[i][j-k]+wk('-',b_base,k) ) {
                        step_right=sw[i][j-k]+wk('-',b_base,k);
                        ki = k;
                    }
                }

                if ( step_down > step_right ) {
                    if ( step_down > step_diag ) {
                        sw[i][j] = Math.max(0,step_down);
                        btrack[i][j] = kd; // positive=vertical
                    }
                    else {
                        sw[i][j] = Math.max(0,step_diag);
                        btrack[i][j] = 0; // 0 = diagonal
                    }
                } else {
                    // step_down < step_right
                    if ( step_right > step_diag ) {
                        sw[i][j] = Math.max(0,step_right);
                        btrack[i][j] = -ki; // negative = horizontal
                    } else {
                        sw[i][j] = Math.max(0,step_diag);
                        btrack[i][j] = 0; // 0 = diagonal
                    }
                }
                sw[i][j] = Math.max(0, Math.max(step_diag,Math.max(step_down,step_right)));
            }
        }

        print(sw,a,b);

        PrimitivePair.Int p = new PrimitivePair.Int();
        double maxscore = 0.0;
        // look for largest score. we use >= combined with the traversal direction
        // to ensure that if two scores are equal, the one closer to diagonal gets picked
        for ( int i = 1 ; i < n+1 ; i++ ) {
            if ( sw[i][m] >= maxscore ) { p.first = i; p.second = m ; maxscore = sw[i][m]; }
        }
        for ( int j = 1 ; j < m+1 ; j++ ) {
            if ( sw[n][j] > maxscore ||
                  sw[n][j] == maxscore && Math.abs(n-j) < Math.abs(p.first-p.second)) {
                p.first = n;
                p.second = j ;
                maxscore = sw[n][j];
            }
        }

        System.out.println("\ni="+p.first+"; j="+p.second);

        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order


        // we will be placing all insertions and deletions into sequence b, so the state are named w/regard
        // to that sequence

        int state = MSTATE;
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

        double [] scores = new double[3];

        List<CigarElement> lce = new ArrayList<CigarElement>(5);

        do {

            int btr = btrack[p.first][p.second];
            int step_left = ( btr < 0 ? -btr : 1);
            int step_up = ( btr > 0 ? btr : 1 );

            // moving left: same base on a, prev base on b = insertion on b:
            scores[ISTATE] = sw[p.first][p.second-step_left] ;
            scores[DSTATE] = sw[p.first - step_up][p.second];
            scores[MSTATE] = sw[p.first-1][p.second-1]; // moving diagonal : match/mismatch

            //            System.out.println("i = " + p.first + " ; j = " + p.second);
            //            System.out.println("s(M)="+scores[MSTATE]+"; s(D)="+scores[DSTATE]+"; s(I)=" + scores[ISTATE]);
            int new_state = findMaxInd(scores,MSTATE);

            int step_length = 1;

            // move to next best location in the sw matrix:
            switch( new_state ) {
                case MSTATE: p.first--; p.second--; break;
                case ISTATE: p.second-=step_left; step_length = step_left; break;
                case DSTATE: p.first-=step_up; step_length = step_up; break;
            }

            // now let's see if the state actually changed:
            if ( new_state == state ) segment_length+=step_length;
            else {
                // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
                CigarOperator o=null;
                switch(state) {
                    case MSTATE: o = CigarOperator.M; break;
                    case ISTATE: o = CigarOperator.I; break;
                    case DSTATE: o = CigarOperator.D; break;
                }
                CigarElement e = new CigarElement(segment_length,o);
                lce.add(e);
                segment_length = step_length;
                state = new_state;
            }
        } while ( scores[state] != 0 );

        // post-process the last segment we are still keeping
        CigarOperator o=null;
        switch(state) {
            case MSTATE: o = CigarOperator.M; break;
            case ISTATE: o = CigarOperator.I; break;
            case DSTATE: o = CigarOperator.D; break;
        }
        CigarElement e = new CigarElement(segment_length,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
        alignment_offset = p.first - p.second;
    }


    private int w(char x, char y) {
        if ( x == y ) return 2; // match
        if ( x == '-' || y == '-' ) return -1; // gap
        return -1; // mismatch
    }

    private double wd ( char x, char y ) {
        if ( x== y ) return 2.0;
        else return -1.0;
    }

    private double wk(char x, char y, int k) {
        return -2.0-k/6; // gap
        // return -1.0 ; // no extension penalty
        // return -1.0-Math.log(k+1); // weak extension penalty
    }

    /** Returns index of the maximum element in array s. If there is a tie, and one of the tied indices is
     * pref_id, then it will be preferred and returned.
     * @param s
     * @param pref_id
     * @return
     */
    private int findMaxInd(int[] s, int pref_id) {
        int imax = 0;
        int maxval = s[0];
        for ( int i = 1; i < s.length ; i++ ) {
            if ( s[i] > maxval || i == pref_id && Math.abs(s[i] - maxval) < 0.0001 ) {
                imax = i;
                maxval = s[i];
            }
        }
        return imax;
    }

    private int findMaxInd(double[] s, int pref_id) {
        int imax = 0;
        double maxval = s[0];
        for ( int i = 1; i < s.length ; i++ ) {
            if ( s[i] > maxval + 0.0001 || i == pref_id && Math.abs(s[i] - maxval) < 0.0001 ) {
                imax = i;
                maxval = s[i];
            }
        }
        return imax;
    }

    private void print(int[][] s) {
        for ( int i = 0 ; i < s.length ; i++) {
            for ( int j = 0; j < s[i].length ; j++ ) {
                System.out.printf(" %4d",s[i][j]);
            }
            System.out.println();
        }
    }

    private void print(double[][] s) {
        for ( int i = 0 ; i < s.length ; i++) {
            for ( int j = 0; j < s[i].length ; j++ ) {
                System.out.printf(" %4g",s[i][j]);
            }
            System.out.println();
        }
    }

    private void print(int[][] s, String a, String b) {

        System.out.print("        ");
        for ( int j = 1 ; j < s[0].length ; j++) System.out.printf(" %4c",b.charAt(j-1)) ;
        System.out.println();

        for ( int i = 0 ; i < s.length ; i++) {
            if ( i > 0 ) System.out.print(a.charAt(i-1));
            else System.out.print(' ');
            System.out.print("  ");
            for ( int j = 0; j < s[i].length ; j++ ) {
                System.out.printf(" %4d",s[i][j]);
            }
            System.out.println();
        }
    }


    private void print(double[][] s, String a, String b) {

        System.out.print("            ");
        for ( int j = 1 ; j < s[0].length ; j++) System.out.printf(" %8c",b.charAt(j-1)) ;
        System.out.println();

        for ( int i = 0 ; i < s.length ; i++) {
            if ( i > 0 ) System.out.print(a.charAt(i-1));
            else System.out.print(' ');
            System.out.print("  ");
            for ( int j = 0; j < s[i].length ; j++ ) {
                System.out.printf(" %8.4f",s[i][j]);
            }
            System.out.println();
        }
    }

    public String toString() {
        StringBuilder b1 = new StringBuilder();
        StringBuilder b2 = new StringBuilder();

        int pos1 = 0;
        int pos2 = 0;
        if ( alignment_offset < 0 ) {
            for (  ; pos2 < -alignment_offset ; pos2++ ) {
                b1.append(' ');
                b2.append(s2.charAt(pos2));
            }
            // now pos2 = -alignment_offset;
        } else {
            for (  ; pos1 < alignment_offset ; pos1++ ) {
                b2.append(' ');
                b1.append(s1.charAt(pos1));
            }
            // now pos1 = alignment_offset
        }

        for ( int i = 0 ; i < getCigar().numCigarElements() ; i++ ) {
            CigarElement ce = getCigar().getCigarElement(i) ;
            switch( ce.getOperator() ) {
                case M:
                    int z = ce.getLength();
                    b1.append(s1, pos1, pos1+z);
                    b2.append(s2, pos2, pos2+z);
                    pos1+=z; pos2+=z;
                    break;
                case I:
                    for ( int k = 0 ; k < ce.getLength() ; k++ ) {
                        b1.append('+');
                        b2.append(s2.charAt(pos2++));
                    }
                    break;
                case D:
                    for ( int k = 0 ; k < ce.getLength() ; k++ ) {
                        b1.append(s1.charAt(pos1++));
                        b2.append('-');
                    }
                    break;
            }
        }
        b1.append(s1,pos1,s1.length());
        b2.append(s2,pos2,s2.length());
        b1.append('\n');
        b1.append(b2);
        b1.append('\n');
        return b1.toString();
    }

    public static void testMe() {
//        String s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
//        String s2 =       "TGTATATAGGGTAAGG";

//        String s1 = "GGTAAGGC";
//        String s2 = "GGTCTCAA";

//        String s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
//        String s2 =       "TGTTAGGGTCTCAAGG";

//        String s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
//        String s2 =             "TAGGGTAAGGCTGATCCATGTACCG" ;

        String s1 =           "ACCTGGTGTATATAGGGTAAGGCTGAT";
        String s2 = "CCGTATCATTACCTGGTGTATATAGG";

//          String s1 = "GGTGTATATAGGGT"  ;
//          String s2 = "TGTTAGGG";
        testMe(s1,s2);
    }

    public static void testMe(String s1, String s2) {

        SWPairwiseAlignment swpa = new SWPairwiseAlignment(s1,s2);


        for ( int i = 0 ; i < swpa.getCigar().numCigarElements() ; i++ ) {
            char c=' ';
            switch ( swpa.getCigar().getCigarElement(i).getOperator() ) {
            case M : c = 'M'; break;
            case D : c = 'D'; break;
            case I : c = 'I'; break;
            }
            System.out.print(c);
            System.out.print(swpa.getCigar().getCigarElement(i).getLength());
        }
        SequencePile sp = new SequencePile(s1);
        sp.addAlignedSequence(s2,false,swpa.getCigar(),swpa.getAlignmentStart2wrt1());
        System.out.println();
        System.out.println(sp.format());

        System.out.println("--------\n"+swpa.toString());        

        //sp.colorprint(false);
    }

    public static void main(String argv[]) {
        if ( argv.length > 0 ) testMe(argv[0],argv[1]);
        else testMe();
    }
}
