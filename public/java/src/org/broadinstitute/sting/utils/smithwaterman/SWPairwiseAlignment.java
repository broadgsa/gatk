/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.smithwaterman;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

import java.util.*;

/**
 * Pairwise discrete smith-waterman alignment
 *
 * ************************************************************************
 * ****                    IMPORTANT NOTE:                             ****
 * ****  This class assumes that all bytes come from UPPERCASED chars! ****
 * ************************************************************************
 *
 * User: asivache
 * Date: Mar 23, 2009
 * Time: 1:54:54 PM
 */
public final class SWPairwiseAlignment {
    private int alignment_offset; // offset of s2 w/respect to s1
    private Cigar alignmentCigar;

    private final Parameters parameters;

    private static final int MSTATE = 0;
    private static final int ISTATE = 1;
    private static final int DSTATE = 2;
    private static final int CLIP = 3;

    protected static boolean cutoff = false;
    private static boolean DO_SOFTCLIP = true;

    /**
     * The SW scoring matrix, stored for debugging purposes if keepScoringMatrix is true
     */
    protected double[] SW = null;

    /**
     * Only for testing purposes in the SWPairwiseAlignmentMain function
     * set to true to keep SW scoring matrix after align call
     */
    protected static boolean keepScoringMatrix = false;

    /**
     * Create a new SW pairwise aligner.
     *
     * @deprecated in favor of constructors using the Parameter or ParameterSet class
     */
    @Deprecated
    public SWPairwiseAlignment(byte[] seq1, byte[] seq2, double match, double mismatch, double open, double extend ) {
        this(seq1, seq2, new Parameters(match, mismatch, open, extend));
    }

    /**
     * Create a new SW pairwise aligner
     *
     * After creating the object the two sequences are aligned with an internal call to align(seq1, seq2)
     *
     * @param seq1 the first sequence we want to align
     * @param seq2 the second sequence we want to align
     * @param parameters the SW parameters to use
     */
    public SWPairwiseAlignment(byte[] seq1, byte[] seq2, Parameters parameters) {
        this.parameters = parameters;
        align(seq1,seq2);
    }

    /**
     * Create a new SW pairwise aligner
     *
     * After creating the object the two sequences are aligned with an internal call to align(seq1, seq2)
     *
     * @param seq1 the first sequence we want to align
     * @param seq2 the second sequence we want to align
     * @param namedParameters the named parameter set to get our parameters from
     */
    public SWPairwiseAlignment(byte[] seq1, byte[] seq2, SWParameterSet namedParameters) {
        this(seq1, seq2, namedParameters.parameters);
    }

    public SWPairwiseAlignment(byte[] seq1, byte[] seq2) {
        this(seq1,seq2,SWParameterSet.ORIGINAL_DEFAULT);
    }

    public Cigar getCigar() { return alignmentCigar ; }

    public int getAlignmentStart2wrt1() { return alignment_offset; }

    public void align(final byte[] a, final byte[] b) {
        final int n = a.length;
        final int m = b.length;
        double [] sw = new double[(n+1)*(m+1)];
        if ( keepScoringMatrix ) SW = sw;
        int [] btrack = new int[(n+1)*(m+1)];

        calculateMatrix(a, b, sw, btrack);
        calculateCigar(n, m, sw, btrack); // length of the segment (continuous matches, insertions or deletions)
    }


    private void calculateMatrix(final byte[] a, final byte[] b, double [] sw, int [] btrack ) {
        final int n = a.length+1;
        final int m = b.length+1;

        //final double MATRIX_MIN_CUTOFF=-1e100;   // never let matrix elements drop below this cutoff
        final double MATRIX_MIN_CUTOFF;   // never let matrix elements drop below this cutoff
        if ( cutoff ) MATRIX_MIN_CUTOFF = 0.0;
        else MATRIX_MIN_CUTOFF = -1e100;

        double [] best_gap_v = new double[m+1];
        Arrays.fill(best_gap_v,-1.0e40);
        int [] gap_size_v = new int[m+1];
        double [] best_gap_h = new double[n+1];
        Arrays.fill(best_gap_h,-1.0e40);
        int [] gap_size_h = new int[n+1];

        // build smith-waterman matrix and keep backtrack info:
        for ( int i = 1, row_offset_1 = 0 ; i < n ; i++ ) { // we do NOT update row_offset_1 here, see comment at the end of this outer loop
            byte a_base = a[i-1]; // letter in a at the current pos

            final int row_offset = row_offset_1 + m;

            // On the entrance into the loop, row_offset_1 is the (linear) offset
            // of the first element of row (i-1) and row_offset is the linear offset of the
            // start of row i

            for ( int j = 1, data_offset_1 = row_offset_1 ; j < m ; j++, data_offset_1++ ) {

                // data_offset_1 is linearized offset of element [i-1][j-1]

                final byte b_base = b[j-1]; // letter in b at the current pos

                // in other words, step_diag = sw[i-1][j-1] + wd(a_base,b_base);
                double step_diag = sw[data_offset_1] + wd(a_base,b_base);

                // optimized "traversal" of all the matrix cells above the current one (i.e. traversing
                // all 'step down' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                // if a gap (length 1) was just opened above, this is the cost of arriving to the current cell:
                double prev_gap = sw[data_offset_1+1]+parameters.w_open;

                best_gap_v[j] += parameters.w_extend; // for the gaps that were already opened earlier, extending them by 1 costs w_extend

                if ( prev_gap > best_gap_v[j] ) {
                    // opening a gap just before the current cell results in better score than extending by one
                    // the best previously opened gap. This will hold for ALL cells below: since any gap
                    // once opened always costs w_extend to extend by another base, we will always get a better score
                    // by arriving to any cell below from the gap we just opened (prev_gap) rather than from the previous best gap
                    best_gap_v[j] = prev_gap;
                    gap_size_v[j] = 1; // remember that the best step-down gap from above has length 1 (we just opened it)
                } else {
                    // previous best gap is still the best, even after extension by another base, so we just record that extension:
                    gap_size_v[j]++;
                }

                final double step_down = best_gap_v[j] ;
                final int kd = gap_size_v[j];

                // optimized "traversal" of all the matrix cells to the left of the current one (i.e. traversing
                // all 'step right' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                final int data_offset = row_offset + j; // linearized offset of element [i][j]
                prev_gap = sw[data_offset-1]+parameters.w_open; // what would it cost us to open length 1 gap just to the left from current cell
                best_gap_h[i] += parameters.w_extend; // previous best gap would cost us that much if extended by another base

                if ( prev_gap > best_gap_h[i] ) {
                    // newly opened gap is better (score-wise) than any previous gap with the same row index i; since
                    // gap penalty is linear with k, this new gap location is going to remain better than any previous ones
                    best_gap_h[i] = prev_gap;
                    gap_size_h[i] = 1;
                } else {
                    gap_size_h[i]++;
                }

                final double step_right = best_gap_h[i];
                final int ki = gap_size_h[i];

                if ( step_down > step_right ) {
                    if ( step_down > step_diag ) {
                        sw[data_offset] = Math.max(MATRIX_MIN_CUTOFF,step_down);
                        btrack[data_offset] = kd ; // positive=vertical
                    } else {
                        sw[data_offset] = Math.max(MATRIX_MIN_CUTOFF,step_diag);
                        btrack[data_offset] = 0; // 0 = diagonal
                    }
                } else {
                    // step_down <= step_right
                    if ( step_right > step_diag ) {
                        sw[data_offset] = Math.max(MATRIX_MIN_CUTOFF,step_right);
                        btrack[data_offset] = -ki; // negative = horizontal
                    } else {
                        sw[data_offset] = Math.max(MATRIX_MIN_CUTOFF,step_diag);
                        btrack[data_offset] = 0; // 0 = diagonal
                    }
                }
            }

            // IMPORTANT, IMPORTANT, IMPORTANT:
            // note that we update this (secondary) outer loop variable here,
            // so that we DO NOT need to update it
            // in the for() statement itself.
            row_offset_1 = row_offset;
        }
    }


    private void calculateCigar(int n, int m, double [] sw, int [] btrack) {
        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order
        int p1 = 0, p2 = 0;

        double maxscore = Double.NEGATIVE_INFINITY; // sw scores are allowed to be negative
        int segment_length = 0; // length of the segment (continuous matches, insertions or deletions)

        // look for largest score. we use >= combined with the traversal direction
        // to ensure that if two scores are equal, the one closer to diagonal gets picked
        for ( int i = 1, data_offset = m+1+m ; i < n+1 ; i++, data_offset += (m+1) ) {
            // data_offset is the offset of [i][m]
            if ( sw[data_offset] >= maxscore ) {
                p1 = i; p2 = m ; maxscore = sw[data_offset];
            }
        }

        for ( int j = 1, data_offset = n*(m+1)+1 ; j < m+1 ; j++, data_offset++ ) {
            // data_offset is the offset of [n][j]
            if ( sw[data_offset] > maxscore || sw[data_offset] == maxscore && Math.abs(n-j) < Math.abs(p1 - p2)) {
                p1 = n;
                p2 = j ;
                maxscore = sw[data_offset];
                segment_length = m - j ; // end of sequence 2 is overhanging; we will just record it as 'M' segment
            }
        }

        List<CigarElement> lce = new ArrayList<CigarElement>(5);

        if ( segment_length > 0 && DO_SOFTCLIP ) {
            lce.add(makeElement(CLIP, segment_length));
            segment_length = 0;
        }

        // we will be placing all insertions and deletions into sequence b, so the states are named w/regard
        // to that sequence

        int state = MSTATE;

        int data_offset = p1*(m+1)+p2;  // offset of element [p1][p2]
        do {
            int btr = btrack[data_offset];

            int new_state;
            int step_length = 1;

            if ( btr > 0 ) {
                new_state = DSTATE;
                step_length = btr;
            } else if ( btr < 0 ) {
                new_state = ISTATE;
                step_length = (-btr);
            } else new_state = MSTATE; // and step_length =1, already set above

            // move to next best location in the sw matrix:
            switch( new_state ) {
                case MSTATE: data_offset -= (m+2); p1--; p2--; break; // move back along the diag in the sw matrix
                case ISTATE: data_offset -= step_length; p2 -= step_length; break; // move left
                case DSTATE: data_offset -= (m+1)*step_length; p1 -= step_length; break; // move up
            }

            // now let's see if the state actually changed:
            if ( new_state == state ) segment_length+=step_length;
            else {
                // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
                lce.add(makeElement(state, segment_length));
                segment_length = step_length;
                state = new_state;
            }
//      next condition is equivalent to  while ( sw[p1][p2] != 0 ) (with modified p1 and/or p2:
        } while ( p1 > 0 && p2 > 0 );

        // post-process the last segment we are still keeping;
        // NOTE: if reads "overhangs" the ref on the left (i.e. if p2>0) we are counting
        // those extra bases sticking out of the ref into the first cigar element if DO_SOFTCLIP is false;
        // otherwise they will be softclipped. For instance,
        // if read length is 5 and alignment starts at offset -2 (i.e. read starts before the ref, and only
        // last 3 bases of the read overlap with/align to the ref), the cigar will be still 5M if
        // DO_SOFTCLIP is false or 2S3M if DO_SOFTCLIP is true.
        // The consumers need to check for the alignment offset and deal with it properly.
        if (DO_SOFTCLIP ) {
            lce.add(makeElement(state, segment_length));
            if ( p2> 0 ) lce.add(makeElement(CLIP, p2));
            alignment_offset = p1 ;
        } else {
            lce.add(makeElement(state, segment_length + p2));
            alignment_offset = p1 - p2;
        }

        Collections.reverse(lce);
        alignmentCigar = AlignmentUtils.consolidateCigar(new Cigar(lce));
    }

    private CigarElement makeElement(int state, int segment_length) {
        CigarOperator o = null;
        switch(state) {
            case MSTATE: o = CigarOperator.M; break;
            case ISTATE: o = CigarOperator.I; break;
            case DSTATE: o = CigarOperator.D; break;
            case CLIP: o = CigarOperator.S; break;
        }
        return new CigarElement(segment_length,o);
    }

    private double wd(byte x, byte y) {
        return (x == y ? parameters.w_match : parameters.w_mismatch);
    }

    public void printAlignment(byte[] ref, byte[] read) {
        printAlignment(ref,read,100);
    }
    
    public void printAlignment(byte[] ref, byte[] read, int width) {
        StringBuilder bread = new StringBuilder();
        StringBuilder bref = new StringBuilder();
        StringBuilder match = new StringBuilder();

        int i = 0;
        int j = 0;

        final int offset = getAlignmentStart2wrt1();

        Cigar cigar = getCigar();

        if ( ! DO_SOFTCLIP ) {

            // we need to go through all the hassle below only if we do not do softclipping;
            // otherwise offset is never negative
            if ( offset < 0 ) {
                for (  ; j < (-offset) ; j++ ) {
                    bread.append((char)read[j]);
                    bref.append(' ');
                    match.append(' ');
                }
                // at negative offsets, our cigar's first element carries overhanging bases
                // that we have just printed above. Tweak the first element to
                // exclude those bases. Here we create a new list of cigar elements, so the original
                // list/original cigar are unchanged (they are unmodifiable anyway!)

                List<CigarElement> tweaked = new ArrayList<CigarElement>();
                tweaked.addAll(cigar.getCigarElements());
                tweaked.set(0,new CigarElement(cigar.getCigarElement(0).getLength()+offset,
                        cigar.getCigarElement(0).getOperator()));
                cigar = new Cigar(tweaked);
            }
        }

        if ( offset > 0 ) { // note: the way this implementation works, cigar will ever start from S *only* if read starts before the ref, i.e. offset = 0
            for (  ; i < getAlignmentStart2wrt1() ; i++ ) {
                bref.append((char)ref[i]);
                bread.append(' ');
                match.append(' ');
            }
        }
        
        for ( CigarElement e : cigar.getCigarElements() ) {
            switch (e.getOperator()) {
                case M :
                    for ( int z = 0 ; z < e.getLength() ; z++, i++, j++  ) {
                        bref.append((i<ref.length)?(char)ref[i]:' ');
                        bread.append((j < read.length)?(char)read[j]:' ');
                        match.append( ( i<ref.length && j < read.length ) ? (ref[i] == read[j] ? '.':'*' ) : ' ' );
                    }
                    break;
                case I :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append('-');
                        bread.append((char)read[j]);
                        match.append('I');
                    }
                    break;
                case S :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append(' ');
                        bread.append((char)read[j]);
                        match.append('S');
                    }
                    break;
                case D:
                    for ( int z = 0 ; z < e.getLength(); z++ , i++ ) {
                        bref.append((char)ref[i]);
                        bread.append('-');
                        match.append('D');
                    }
                    break;
                default:
                    throw new StingException("Unexpected Cigar element:" + e.getOperator());
            }
        }
        for ( ; i < ref.length; i++ ) bref.append((char)ref[i]);
        for ( ; j < read.length; j++ ) bread.append((char)read[j]);

        int pos = 0 ;
        int maxlength = Math.max(match.length(),Math.max(bread.length(),bref.length()));
        while ( pos < maxlength ) {
            print_cautiously(match,pos,width);
            print_cautiously(bread,pos,width);
            print_cautiously(bref,pos,width);
            System.out.println();
            pos += width;
        }
    }

    /** String builder's substring is extremely stupid: instead of trimming and/or returning an empty
     * string when one end/both ends of the interval are out of range, it crashes with an
     * exception. This utility function simply prints the substring if the interval is within the index range
     * or trims accordingly if it is not.
     * @param s
     * @param start
     * @param width
     */
    private static void print_cautiously(StringBuilder s, int start, int width) {
        if ( start >= s.length() ) {
            System.out.println();
            return;
        }
        int end = Math.min(start+width,s.length());
        System.out.println(s.substring(start,end));
    }
}
