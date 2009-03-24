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
        int segment_length = 1; // length of the segment (continuous matches, insertions or deletions)

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
        CigarElement e = new CigarElement(segment_length-1,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
        alignment_offset = p.first;
    }

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
                for ( int k = 1 ; k < j ; k++ ) step_down = Math.max(step_down,sw[i][j-k]+wk('-',b_base,k));

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
        int segment_length = 1; // length of the segment (continuous matches, insertions or deletions)

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
        CigarElement e = new CigarElement(segment_length-1,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
        alignment_offset = p.first;
    }


    private int w(char x, char y) {
        if ( x == y ) return 2; // match
        if ( x == '-' || y == '-' ) return -1; // gap
        return -1; // mismatch
    }

    private double wd ( char x, char y ) {
        if ( x== y ) return 2.0;
        else return -1;
    }

    private double wk(char x, char y, int k) {
        return -2.0-k; // gap
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

    public static void testMe() {
//        String s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
//        String s2 =       "TGTATATAGGGTAAGG";

        String s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
        String s2 =       "TGTTAGGGTCTCAAGG";

        testMe(s1,s2);
    }

    public static void testMe(String s1, String s2) {

        SWPairwiseAlignment swpa = new SWPairwiseAlignment(s1,s2);

        SequencePile sp = new SequencePile(s1);
        sp.addAlignedSequence(s2,false,swpa.getCigar(),swpa.getAlignmentStart2wrt1());

        for ( int i = 0 ; i < swpa.getCigar().numCigarElements() ; i++ ) {
            System.out.print(swpa.getCigar().getCigarElement(i).getLength());
            char c=' ';
            switch ( swpa.getCigar().getCigarElement(i).getOperator() ) {
            case M : c = 'M'; break;
            case D : c = 'D'; break;
            case I : c = 'I'; break;
            }
            System.out.print(c);
        }
        System.out.println();
        System.out.println(sp.format());

        //sp.colorprint(false);
    }

    public static void main(String argv[]) {
        if ( argv.length > 0 ) testMe(argv[0],argv[1]);
        else testMe();
    }
}
