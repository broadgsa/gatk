/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.Cigar;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

import org.broadinstitute.sting.utils.collections.PrimitivePair;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 23, 2009
 * Time: 1:54:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class SWPairwiseAlignment {
    private int alignment_offset; // offset of s2 w/respect to s1
    private Cigar alignmentCigar;

    private final double w_match;
    private final double w_mismatch;
    private final double w_open;
    private final double w_extend;

    private static final int MSTATE = 0;
    private static final int ISTATE = 1;
    private static final int DSTATE = 2;

    // ************************************************************************
    // ****                    IMPORTANT NOTE:                             ****
    // ****  This class assumes that all bytes come from UPPERCASED chars! ****
    // ************************************************************************
    public SWPairwiseAlignment(byte[] seq1, byte[] seq2, double match, double mismatch, double open, double extend ) {
        w_match = match;
        w_mismatch = mismatch;
        w_open = open;
        w_extend = extend;
        align(seq1,seq2);
    }

    public SWPairwiseAlignment(byte[] seq1, byte[] seq2) {
        this(seq1,seq2,1.0,-1.0/3.0,-1.0-1.0/3.0,-1.0/3.0); // match=1, mismatch = -1/3, gap=-(1+k/3)
    }

    public Cigar getCigar() { return alignmentCigar ; }

    public int getAlignmentStart2wrt1() { return alignment_offset; }

    public void align(final byte[] a, final byte[] b) {
            final int n = a.length;
            final int m = b.length;
            double [][] sw = new double[n+1][m+1];

            int [][] btrack = new int[n+1][m+1];

            // build smith-waterman matrix and keep backtrack info:
            for ( int i = 1 ; i < n+1 ; i++ ) {
                byte a_base = a[i-1]; // letter in a at the current pos
                for ( int j = 1 ; j < m+1 ; j++ ) {
                    final byte b_base = b[j-1]; // letter in b at the current pos
                    double step_diag = sw[i-1][j-1] + wd(a_base,b_base);
                    double step_down = 0.0 ;
                    int kd = 0;
                    for ( int k = 1 ; k < i ; k++ ) {
                        final double gap_penalty = wk(k);
                        if ( step_down < sw[i-k][j]+gap_penalty ) {
                            step_down=sw[i-k][j]+gap_penalty;
                            kd = k;
                        }
                    }

                    double step_right = 0;
                    int ki = 0;
                    for ( int k = 1 ; k < j ; k++ ) {
                        final double gap_penalty = wk(k);
                        if ( step_right < sw[i][j-k]+gap_penalty ) {
                            step_right=sw[i][j-k]+gap_penalty;
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

            List<CigarElement> lce = new ArrayList<CigarElement>(5);

            do {

                int btr = btrack[p.first][p.second];
                int step_left = ( btr < 0 ? -btr : 1);
                int step_up = ( btr > 0 ? btr : 1 );

                int new_state;
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

    private double wd(byte x, byte y) {
        return (x == y ? w_match : w_mismatch);
    }

    private double wk(int k) {
        return w_open+(k-1)*w_extend; // gap
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

}
