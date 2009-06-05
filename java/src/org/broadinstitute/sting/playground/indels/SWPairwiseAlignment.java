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

    private double w_match;
    private double w_mismatch;
    private double w_open;
    private double w_extend;

    private int best_mm; // mismatch count
    private static final int IMPOSSIBLE = 1000000000;
    private static final int MSTATE = 0;
    private static final int ISTATE = 1;
    private static final int DSTATE = 2;

    public SWPairwiseAlignment(String seq1, String seq2, int id1, int id2, double match, double mismatch, double open, double extend ) {
        s1 = seq1;
        s2 = seq2;
        i1 = id1;
        i2 = id2;
        w_match = match;
        w_mismatch = mismatch;
        w_open = open;
        w_extend = extend;
        best_mm = IMPOSSIBLE;
        //next_mm = IMPOSSIBLE;
        align4(s1,s2);
    }

    public SWPairwiseAlignment(String seq1, String seq2, int id1, int id2) {
        this(seq1,seq2,id1,id2,1.0,-1.0/3.0,-1.0-1.0/3.0,-1.0/3.0); // match=1, mismatch = -1/3, gap=-(1+k/3)
    }

    /** Initializes the alignment with pair of sequences (that will be immediately aligned) and
     * sets their external ids to -1. Such un-annotated pairwise alignment can not be added to MultipleAlignment.
     *
     */
    public SWPairwiseAlignment(String seq1, String seq2) {
        this(seq1,seq2,-1,-1);
    }

    public SWPairwiseAlignment(String seq1, String seq2, double match, double mismatch, double open, double extend) {
        this(seq1,seq2,-1,-1,match,mismatch,open, extend);
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

//        print(sw,a,b);

        PrimitivePair.Int p = new PrimitivePair.Int();
        int maxscore = 0;
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

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
                segment_length = m - j; // end of sequence 2 is overhanging; we will just record it as 'M' segment
            }
        }

//        System.out.println("\ni="+p.first+"; j="+p.second);

        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order


        // we will be placing all insertions and deletions into sequence b, so the state are named w/regard
        // to that sequence

        int state = MSTATE;

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

        alignment_offset = p.first - p.second;
        segment_length+=p.second;
        CigarElement e = new CigarElement(segment_length,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
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

//        print(sw,a,b);

        PrimitivePair.Int p = new PrimitivePair.Int();
        double maxscore = 0.0;
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

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
                segment_length = m - j; // end of sequence 2 is overhanging; we will just record it as 'M' segment
            }
        }

//        System.out.println("\ni="+p.first+"; j="+p.second);

        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order


        // we will be placing all insertions and deletions into sequence b, so the state are named w/regard
        // to that sequence

        int state = MSTATE;

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
        alignment_offset = p.first - p.second;
        segment_length+=p.second;
        CigarElement e = new CigarElement(segment_length,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
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

//        print(sw,a,b);

        PrimitivePair.Int p = new PrimitivePair.Int();
        double maxscore = 0.0;
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

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
                segment_length = m - j ; // end of sequence 2 is overhanging; we will just record it as 'M' segment
            }
        }

        System.out.println("\ni="+p.first+"; j="+p.second);

        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order


        // we will be placing all insertions and deletions into sequence b, so the state are named w/regard
        // to that sequence

        int state = MSTATE;

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
        alignment_offset = p.first - p.second;
        segment_length+=p.second;
        CigarElement e = new CigarElement(segment_length,o);
        lce.add(e);
        Collections.reverse(lce);
        alignmentCigar = new Cigar(lce);
    }

    public void align4(String a, String b) {
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

//            print(sw,a,b);

            PrimitivePair.Int p = new PrimitivePair.Int();
            double maxscore = 0.0;
            int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

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
                    segment_length = m - j ; // end of sequence 2 is overhanging; we will just record it as 'M' segment
                }
            }

 //           System.out.println("\ni="+p.first+"; j="+p.second);

            // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order


            // we will be placing all insertions and deletions into sequence b, so the state are named w/regard
            // to that sequence

            int state = MSTATE;

            double [] scores = new double[3];

            List<CigarElement> lce = new ArrayList<CigarElement>(5);

            do {

                int btr = btrack[p.first][p.second];
                int step_left = ( btr < 0 ? -btr : 1);
                int step_up = ( btr > 0 ? btr : 1 );

                int new_state = -1;
                if ( btr > 0 ) new_state = DSTATE;
                else if ( btr < 0 ) new_state = ISTATE;
                else new_state = MSTATE;

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
            } while ( sw[p.first][p.second] != 0 );

            // post-process the last segment we are still keeping
            CigarOperator o=null;
            switch(state) {
                case MSTATE: o = CigarOperator.M; break;
                case ISTATE: o = CigarOperator.I; break;
                case DSTATE: o = CigarOperator.D; break;
            }
            alignment_offset = p.first - p.second;
            segment_length+=p.second;
            CigarElement e = new CigarElement(segment_length,o);
            lce.add(e);
            Collections.reverse(lce);
            alignmentCigar = new Cigar(lce);
    }

    private int w(char x, char y) {
        if ( x == y ) return 2; // match
        if ( x == '-' || y == '-' ) return -1; // gap
        return -1; // mismatch
    }

    private double wd ( char x, char y ) {
        if ( x== y ) return w_match;
        else return w_mismatch;
    }

    private double wk(char x, char y, int k) {
        return w_open+(k-1)*w_extend; // gap
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

        System.out.print("");
        for ( int j = 1 ; j < s[0].length ; j++) System.out.printf(" %4c",b.charAt(j-1)) ;
        System.out.println();

        for ( int i = 0 ; i < s.length ; i++) {
            if ( i > 0 ) System.out.print(a.charAt(i-1));
            else System.out.print(' ');
            System.out.print("  ");
            for ( int j = 0; j < s[i].length ; j++ ) {
                System.out.printf(" %2.1f",s[i][j]);
            }
            System.out.println();
        }
    }

    public String toString() {
        StringBuilder bmm = new StringBuilder();
        StringBuilder b1 = new StringBuilder();
        StringBuilder b2 = new StringBuilder();

        int pos1 = 0;
        int pos2 = 0;
        if ( alignment_offset < 0 ) {
            for (  ; pos2 < -alignment_offset ; pos2++ ) {
                b1.append(' ');
                b2.append(s2.charAt(pos2));
                bmm.append(' ');
            }
            // now pos2 = -alignment_offset;
        } else {
            for (  ; pos1 < alignment_offset ; pos1++ ) {
                b2.append(' ');
                b1.append(s1.charAt(pos1));
                bmm.append(' ');
            }
            // now pos1 = alignment_offset
        }
/* debug prints: */
//        System.out.println(AlignmentUtils.toString(getCigar()));
//        System.out.println("seq1l="+s1.length()+"; seq2l=" + s2.length());
//        System.out.println("offset="+alignment_offset);
//        System.out.println("pos1="+pos1+"; pos2=" + pos2);
/**/
        for ( int i = 0 ; i < getCigar().numCigarElements() ; i++ ) {
            CigarElement ce = getCigar().getCigarElement(i) ;
            switch( ce.getOperator() ) {
                case M:
                    int z = ( i == 0 ? pos2 : 0); // if we are in the first element and seq overhangs to the left,
                                                  // start inside the first segment, at the position where actual matches begin
                    // check separately for pos1 < s1.length() since seq2 is allowed to overhang beyond seq1's end
                    for (  ; z < ce.getLength() && pos1 < s1.length() ; z++ ) {
//                        System.out.println("pos1="+pos1+"; pos2="+pos2+"; k="+z);
                        if ( Character.toUpperCase(s1.charAt(pos1)) !=
                                Character.toUpperCase(s2.charAt(pos2)) ) bmm.append('*');
                        else bmm.append(' ');
                        b1.append(s1.charAt(pos1++));
                        b2.append(s2.charAt(pos2++));
                    }
                    break;
                case I:
                    for ( int k = 0 ; k < ce.getLength() ; k++ ) {
                        b1.append('+');
                        bmm.append('+');
                        b2.append(s2.charAt(pos2++));
                    }
                    break;
                case D:
                    for ( int k = 0 ; k < ce.getLength() ; k++ ) {
                        b1.append(s1.charAt(pos1++));
                        bmm.append('-');
                        b2.append('-');
                    }
                    break;
            }
        }

        bmm.append('\n');
        b1.append(s1,pos1,s1.length());
        bmm.append(b1);
        bmm.append('\n');
        b2.append(s2,pos2,s2.length());
        bmm.append(b2);
        bmm.append('\n');

        return bmm.toString();
    }

    public static void testMe() {
 /*       String s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
        String s2 =       "TGTATATAGGGTAAGG";

        testMe(s1,s2);

        s1 = "GGTAAGGC";
        s2 = "GGTCTCAA";

        testMe(s1,s2);

        s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
        s2 =       "TGTTAGGGTCTCAAGG";
        testMe(s1,s2);


        s1 = "ACCTGGTGTATATAGGGTAAGGCTGAT";
        s2 =             "TAGGGTAAGGCTGATCCATGTACCG" ;
        testMe(s1,s2);

        s1 =           "ACCTGGTGTATATAGGGTAAGGCTGAT";
        s2 = "CCGTATCATTACCTGGTGTATATAGG";
        testMe(s1,s2);

        s1 = "GGTGTATATAGGGT"  ;
        s2 = "TGTTAGGG";
        testMe(s1,s2);

        s1 = "AGACAGAGAGAAGG";
        s2 = "AGACAGAGAAGG";
        testMe(s1,s2);
   */
 //       String s1 = "CCAGCACACAGGTATCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTGTTTTTTGA";
 //       String s2 = "CCAGCACACATCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTGTTTTTTGA";

//        String s1 = "CCCATCTGTCTCCAATCTGCTGTTTTCCAAAAATTAGGGAACTTCAGTTTTCCCTTTGATACTCTGTATTTCTACCAACCACAACGCCAGGGCTGTCCTGCTTCTACAAGTGACAATGACAAATATAGGCCTGAAGGAAGATG";
//        String s2 = "AAAATCTGTTTCCAATCTACTGTTTTCCAAAAATTAGGGAAGTTCAGTTTTCCCTTTGATACTCTGTTTCTACCAATCC";
        String s1 = "CCCATCTGTCTCCAATCTGCTGTTTTCCAAAAATTAGGGAACTTCAGTTTTCCCTTTGATACTCTGTATTTCTACCAACCACAACGCCAGGGCTGTCCTGCTTCTACAAGTGACAATGACAAATATAGGCCTGAAGGAAGATG";
        String s2 = "AAAATCTGTCTCCAATCTACTGTTTTCCAAAAATTAGGGAAGTTCAGTTTTCCCTTTGATACTCTGTTTCTACCAATCC";
           testMe(s1,s2);
    }

    public static void testMe(String s1, String s2) {

        SWPairwiseAlignment swpa = new SWPairwiseAlignment(s1,s2,3.0,-1.0,-4,-0.5);

        System.out.println(AlignmentUtils.toString(swpa.getCigar()));
//       SequencePile sp = new SequencePile(s1);
//        sp.addAlignedSequence(s2,false,swpa.getCigar(),swpa.getAlignmentStart2wrt1());
//        System.out.println();
//        System.out.println(sp.format());

        System.out.println("--------\n"+swpa.toString());        

        //sp.colorprint(false);
    }

    public static void main(String argv[]) {
        if ( argv.length > 0 ) testMe(argv[0],argv[1]);
        else testMe();
    }
}
