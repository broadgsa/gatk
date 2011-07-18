/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.*;

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
    private static final int CLIP = 3;

    private static boolean cutoff = false;
    private static boolean DO_SOFTCLIP = true;

    double[] SW;

//    private double [] best_gap_v ;
//    private int [] gap_size_v ;
//    private double [] best_gap_h ;
//    private int [] gap_size_h ;


 //   private static double [][] sw = new double[500][500];
 //   private static int [][] btrack = new int[500][500];

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
        double [] sw = new double[(n+1)*(m+1)];
        SW = sw;
        int [] btrack = new int[(n+1)*(m+1)];

//        best_gap_v = new double[m+1];
//        Arrays.fill(best_gap_v,-1.0e40);
//        gap_size_v = new int[m+1];
//        best_gap_h = new double[n+1];
//        Arrays.fill(best_gap_h,-1.0e40);
//        gap_size_h = new int[n+1];

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
                double prev_gap = sw[data_offset_1+1]+w_open;

                best_gap_v[j] += w_extend; // for the gaps that were already opened earlier, extending them by 1 costs w_extend

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

/*
                for ( int k = 1, data_offset_k = data_offset_1+1 ; k < i ; k++, data_offset_k -= m ) {
                    // data_offset_k is linearized offset of element [i-k][j]
                    // in other words, trial = sw[i-k][j]+gap_penalty:
                    final double trial = sw[data_offset_k]+wk(k);
                    if ( step_down < trial ) {
                        step_down=trial;
                        kd = k;
                    }
                }
*/

                // optimized "traversal" of all the matrix cells to the left of the current one (i.e. traversing
                // all 'step right' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                final int data_offset = row_offset + j; // linearized offset of element [i][j]
                prev_gap = sw[data_offset-1]+w_open; // what would it cost us to open length 1 gap just to the left from current cell
                best_gap_h[i] += w_extend; // previous best gap would cost us that much if extended by another base

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

/*
                for ( int k = 1, data_offset = row_offset+j-1 ; k < j ; k++, data_offset-- ) {
                    // data_offset is linearized offset of element [i][j-k]
                    // in other words, step_right=sw[i][j-k]+gap_penalty;
                    final double trial = sw[data_offset]+wk(k);
                    if ( step_right < trial ) {
                        step_right=trial;
                        ki = k;
                    }
                }

                final int data_offset = row_offset + j; // linearized offset of element [i][j]
*/


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

//                sw[data_offset] = Math.max(0, Math.max(step_diag,Math.max(step_down,step_right)));
            }

            // IMPORTANT, IMPORTANT, IMPORTANT:
            // note that we update this (secondary) outer loop variable here,
            // so that we DO NOT need to update it
            // in the for() statement itself.
            row_offset_1 = row_offset;
        }
//            print(sw,a,b);
    }


    private void calculateCigar(int n, int m, double [] sw, int [] btrack) {
        // p holds the position we start backtracking from; we will be assembling a cigar in the backwards order
        //PrimitivePair.Int p = new PrimitivePair.Int();
        int p1 = 0, p2 = 0;

        double maxscore = 0.0;
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
//                maxscore = sw[n][j];
                maxscore = sw[data_offset];
                segment_length = m - j ; // end of sequence 2 is overhanging; we will just record it as 'M' segment
            }
        }
//        System.out.println("  Found max score="+maxscore+" at p1="+p1+ " p2="+p2);

        List<CigarElement> lce = new ArrayList<CigarElement>(5);

        if ( segment_length > 0 && DO_SOFTCLIP ) {
            lce.add(makeElement(CLIP, segment_length));
            segment_length = 0;
        }

        // we will be placing all insertions and deletions into sequence b, so the states are named w/regard
        // to that sequence

        int state = MSTATE;

        int data_offset = p1*(m+1)+p2;  // offset of element [p1][p2]
 //       System.out.println("Backtracking: starts at "+p1+":"+p2+" ("+sw[data_offset]+")");
        do {
//            int btr = btrack[p1][p2];
            int btr = btrack[data_offset];

            int new_state;
            int step_length = 1;

 //           System.out.print(" backtrack value: "+btr);

            if ( btr > 0 ) {
                new_state = DSTATE;
                step_length = btr;
            } else if ( btr < 0 ) {
                new_state = ISTATE;
                step_length = (-btr);
            } else new_state = MSTATE; // and step_length =1, already set above


            // move to next best location in the sw matrix:
            switch( new_state ) {
                case MSTATE: data_offset -= (m+2); p1--; p2--; break; // move back along the diag in th esw matrix
                case ISTATE: data_offset -= step_length; p2 -= step_length; break; // move left
                case DSTATE: data_offset -= (m+1)*step_length; p1 -= step_length; break; // move up
            }
  //          System.out.println("; backtracked to p1="+p1+" p2="+p2);
  /*
            switch( new_state ) {
                case MSTATE: System.out.println("  diag (match) to "+ sw[data_offset]); break; // equivalent to p1--; p2--
                case ISTATE: System.out.println("  left (insertion, "+step_length+") to "+ sw[data_offset]); break; // equivalent to p2-=step_length;
                case DSTATE: System.out.println("    up (deletion, "+step_length+") to "+ sw[data_offset]); break; // equivalent to p1 -= step_up
            }
   */
            // now let's see if the state actually changed:
            if ( new_state == state ) segment_length+=step_length;
            else {
//                System.out.println(" emitting "+segment_length+makeElement(state,segment_length).getOperator().toString());
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
        alignmentCigar = new Cigar(lce);

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

    private void print(double[] s, byte[] a, byte[] b) {
           int n = a.length+1;
           int m = b.length+1;
           System.out.print("         ");
           for ( int j = 1 ; j < m ; j++) System.out.printf(" %5c",(char)b[j-1]) ;
           System.out.println();

           for ( int i = 0, row_offset = 0 ; i < n ; i++, row_offset+=m) {
               if ( i > 0 ) System.out.print((char)a[i-1]);
               else System.out.print(' ');
               System.out.print("  ");
               for ( int j = 0; j < m ; j++ ) {
                   System.out.printf(" %5.1f",s[row_offset+j]);
               }
               System.out.println();
           }
       }

    static void printAlignment(SWPairwiseAlignment a, byte[] ref, byte[] read) {
        printAlignment(a,ref,read,100);
    }
    
    static void printAlignment(SWPairwiseAlignment a, byte[] ref, byte[] read, int width) {
        StringBuilder bread = new StringBuilder();
        StringBuilder bref = new StringBuilder();
        StringBuilder match = new StringBuilder();

        int i = 0;
        int j = 0;

        final int offset = a.getAlignmentStart2wrt1();

        Cigar cigar = a.getCigar();

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
            for (  ; i < a.getAlignmentStart2wrt1() ; i++ ) {
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

//    BELOW: main() method for testing; old implementations of the core methods are commented out below;
//           uncomment everything through the end of the file if benchmarking of new vs old implementations is needed.

    public static void main(String argv[]) {
//        String ref="CACGAGCATATGTGTACATGAATTTGTATTGCACATGTGTTTAATGCGAACACGTGTCATGTGTATGTGTTCACATGCATGTGTGTCT";
//        String read =   "GCATATGTTTACATGAATTTGTATTGCACATGTGTTTAATGCGAACACGTGTCATGTGTGTGTTCACATGCATGTG";

        String ref = null;
        String read = null;

        Map<String,List<String>> args = processArgs(argv);

        List<String> l = args.get("SEQ");
        args.remove("SEQ");
        if ( l == null ) {
            System.err.println("SEQ argument is missing. Two input sequences must be provided");
            System.exit(1);
        }
        if ( l.size() != 2 ) {
            System.err.println("Two input sequences (SEQ arguments) must be provided. Found "+l.size()+" instead");
            System.exit(1);
        }

        ref = l.get(0);
        read = l.get(1);

        Double m = extractSingleDoubleArg("MATCH",args);
        Double mm = extractSingleDoubleArg("MISMATCH",args);
        Double open = extractSingleDoubleArg("OPEN",args);
        Double ext = extractSingleDoubleArg("EXTEND",args);

        Boolean reverse = extractSingleBooleanArg("REVERSE",args);
        if ( reverse != null && reverse.booleanValue() == true ) {
            ref = Utils.reverse(ref);
            read = Utils.reverse(read);
        }

        Boolean print_mat = extractSingleBooleanArg("PRINT_MATRIX",args);
        Boolean cut = extractSingleBooleanArg("CUTOFF",args);
        if ( cut != null ) SWPairwiseAlignment.cutoff = cut;

        if ( args.size() != 0 ) {
            System.err.println("Unknown argument on the command line: "+args.keySet().iterator().next());
            System.exit(1);
        }

        double w_match;
        double w_mismatch;
        double w_open;
        double w_extend;

        w_match = (m == null ? 30.0 : m.doubleValue());
        w_mismatch = (mm == null ? -10.0 : mm.doubleValue());
        w_open = (open == null ? -10.0 : open.doubleValue());
        w_extend = (ext == null ? -2.0 : ext.doubleValue());


        SWPairwiseAlignment a = new SWPairwiseAlignment(ref.getBytes(),read.getBytes(),w_match,w_mismatch,w_open,w_extend);

        System.out.println("start="+a.getAlignmentStart2wrt1()+", cigar="+a.getCigar()+
                " length1="+ref.length()+" length2="+read.length());


        System.out.println();
        printAlignment(a,ref.getBytes(),read.getBytes());

        System.out.println();
        if ( print_mat != null && print_mat == true ) {
            a.print(a.SW,ref.getBytes(),read.getBytes());
        }
    }


    static Pair<String,Integer> getArg(String prefix, String argv[], int i) {
        String arg = null;
        if ( argv[i].startsWith(prefix) ) {
            arg = argv[i].substring(prefix.length());
            if( arg.length() == 0 ) {
                i++;
                if ( i < argv.length ) arg = argv[i];
                else {
                    System.err.println("No value found after " + prefix + " argument tag");
                    System.exit(1);
                }
            }
            i++;
        }
        return new Pair<String,Integer>(arg,i);
    }

    static Map<String,List<String>> processArgs(String argv[]) {
        Map<String,List<String>> args = new HashMap<String,List<String>>();

        for ( int i = 0; i < argv.length ; i++ ) {
            String arg = argv[i];
            int pos = arg.indexOf('=');
            if ( pos < 0 ) {
                System.err.println("Argument "+arg+" is not of the form <ARG>=<VAL>");
                System.exit(1);
            }
            String val = arg.substring(pos+1);
            if ( val.length() == 0 ) {
                // there was a space between '=' and the value
                i++;
                if ( i < argv.length ) val = argv[i];
                else {
                    System.err.println("No value found after " + arg + " argument tag");
                    System.exit(1);
                }
            }
            arg = arg.substring(0,pos);

            List<String> l = args.get(arg);
            if ( l == null ) {
                l = new ArrayList<String>();
                args.put(arg,l);
            }
            l.add(val);
        }
        return args;
    }

    static Double extractSingleDoubleArg(String argname, Map<String,List<String>> args) {
        List<String> l = args.get(argname);
        args.remove(argname);
        if ( l == null ) return null;

        if ( l.size() > 1 ) {
            System.err.println("Only one "+argname+" argument is allowed");
            System.exit(1);
        }
        double d=0;
        try {
            d = Double.parseDouble(l.get(0));
        } catch ( NumberFormatException e) {
            System.err.println("Can not parse value provided for "+argname+" argument ("+l.get(0)+")");
            System.exit(1);
        }
        System.out.println("Argument "+argname+" set to "+d);
        return new Double(d);
    }


    static Boolean extractSingleBooleanArg(String argname, Map<String,List<String>> args) {
        List<String> l = args.get(argname);
        args.remove(argname);
        if ( l == null ) return null;

        if ( l.size() > 1 ) {
            System.err.println("Only one "+argname+" argument is allowed");
            System.exit(1);
        }
        if ( l.get(0).equals("true") ) return new Boolean(true);
        if ( l.get(0).equals("false") ) return new Boolean(false);
        System.err.println("Can not parse value provided for "+argname+" argument ("+l.get(0)+"); true/false are allowed");
        System.exit(1);
        return null;
    }

/* ##############################################
    public SWPairwiseAlignment(byte[] seq1, byte[] seq2, double match, double mismatch, double open, double extend, boolean runOld ) {
        w_match = match;
        w_mismatch = mismatch;
        w_open = open;
        w_extend = extend;
        if ( runOld ) align_old(seq1,seq2);
        else align(seq1,seq2);
    }

    public SWPairwiseAlignment(byte[] seq1, byte[] seq2, boolean runOld) {
        this(seq1,seq2,1.0,-1.0/3.0,-1.0-1.0/3.0,-1.0/3.0,runOld); // match=1, mismatch = -1/3, gap=-(1+k/3)
    }

    public void align_old(final byte[] a, final byte[] b) {
        final int n = a.length;
        final int m = b.length;
        double [] sw = new double[(n+1)*(m+1)];
        int [] btrack = new int[(n+1)*(m+1)];
        calculateMatrix_old(a, b, sw, btrack);
        calculateCigar(n, m, sw, btrack); // length of the segment (continuous matches, insertions or deletions)
    }

    private void calculateMatrix_old(final byte[] a, final byte[] b, double [] sw, int [] btrack ) {
        final int n = a.length+1;
        final int m = b.length+1;

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
                int kd = 0;

                double step_down = 0;

                for ( int k = 1, data_offset_k = data_offset_1+1 ; k < i ; k++, data_offset_k -= m ) {
                    // data_offset_k is linearized offset of element [i-k][j]
                    // in other words, trial = sw[i-k][j]+gap_penalty:
                    final double trial = sw[data_offset_k]+wk(k);
                    if ( step_down < trial ) {
                        step_down=trial;
                        kd = k;
                    }
                }

                int ki = 0;

                // optimized "traversal" of all the matrix cells to the left of the current one (i.e. traversing
                // all 'step right' events that would end in the current cell. The optimized code
                // does exactly the same thing as the commented out loop below. IMPORTANT:
                // the optimization works ONLY for linear w(k)=wopen+(k-1)*wextend!!!!

                double step_right = 0;

                for ( int k = 1, data_offset = row_offset+j-1 ; k < j ; k++, data_offset-- ) {
                    // data_offset is linearized offset of element [i][j-k]
                    // in other words, step_right=sw[i][j-k]+gap_penalty;
                    final double trial = sw[data_offset]+wk(k);
                    if ( step_right < trial ) {
                        step_right=trial;
                        ki = k;
                    }
                }

                final int data_offset = row_offset + j; // linearized offset of element [i][j]

                if ( step_down > step_right ) {
                    if ( step_down > step_diag ) {
                        sw[data_offset] = Math.max(0,step_down);
                        btrack[data_offset] = kd ; // positive=vertical
                    } else {
                        sw[data_offset] = Math.max(0,step_diag);
                        btrack[data_offset] = 0; // 0 = diagonal
                    }
                } else {
                    // step_down <= step_right
                    if ( step_right > step_diag ) {
                        sw[data_offset] = Math.max(0,step_right);
                        btrack[data_offset] = -ki; // negative = horizontal
                    } else {
                        sw[data_offset] = Math.max(0,step_diag);
                        btrack[data_offset] = 0; // 0 = diagonal
                    }
                }

//                sw[data_offset] = Math.max(0, Math.max(step_diag,Math.max(step_down,step_right)));
            }

            // IMPORTANT, IMPORTANT, IMPORTANT:
            // note that we update this (secondary) outer loop variable here,
            // so that we DO NOT need to update it
            // in the for() statement itself.
            row_offset_1 = row_offset;
        }
//            print(sw,a,b);
    }
#####################
END COMMENTED OUT SECTION
*/

}
