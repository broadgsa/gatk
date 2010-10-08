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

package org.broadinstitute.sting.playground.tools;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.samtools.*;

import java.io.File;
import java.util.List;
import java.util.ArrayList;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.playground.utils.ParallelSAMIterator;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Aug 27, 2009
 * Time: 3:53:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class PairMaker extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Reconstitutes mate pairs from alignments"+
                " for individual fragment ends. Individual end alignment files are expected to be sorted by"+
            " read name. Multiple alignments are allowed and in this case the single best pairing will be selected.";
    @Option(shortName="I1",
            doc="Input file (bam or sam) with alignments for end 1 (first mate in a pair).",
            optional=false)
    public File IN1 = null;
    @Option(shortName="I2",
            doc="Input file (bam or sam) with alignments for end 2 (second mate in a pair).",
            optional=false)
    public File IN2 = null;
    @Option(shortName="O",optional=false, doc="Output file to write found/selected pairs into.")
        public File OUTPUT = null;
    @Option(shortName="F",doc="Turns on a 'filter' mode: only records/pairs passing the filter will be written "+
                "into the output file. Filter condition is a logical combination (using parentheses for grouping, "+
                "& for AND, | for OR, = for specifying values) of the primitives listed below. Primitives that end with "+
                " 1 or 2 apply specifically to pair end 1 or 2, respectively; when, in addition to primitives <P>1 and <P>2, " +
                " a primitive <P> is also defined, it is always interpreted as <P>1 & <P>2. Primitives: PAIRSONLY "+
                "(print only pairs=both ends mapped), MINQ1=<N>, MINQ2=<N>, MINQ=<N> (minimum alignment quality if read is mapped) "+
                " for MINQ1=<N> & MINQ2=<N>), ERRDIFF1=<N>, ERRDIFF2=<N> " +
                "(when next-best alignment(s) are available",
            optional = true)
        public String FILTER = null;
    @Option(shortName="Q", optional=true, doc="Minimum mapping quality required on both ends in order to accept the pair.")
        public Integer MINQ = 20;

    public static int INFINITY = 1000000000;

    // we will collect some stats along the way:
    private int fragments_seen = 0;
    private int end1_missing = 0;
    private int end2_missing = 0;
    private int end1_unmapped = 0;
    private int end2_unmapped = 0;
    private int both_unmapped = 0;
    private int both_mapped = 0;
    private int both_unique = 0;
    private int proper_pair = 0;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new PairMaker().instanceMain(argv));
    }

    protected int doWork() {

        SAMFileReader end1Reader = new SAMFileReader(IN1);
        SAMFileReader end2Reader = new SAMFileReader(IN2);

        SAMFileHeader h = checkHeaders(end1Reader.getFileHeader(), end2Reader.getFileHeader() );

        SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, OUTPUT);

        ParallelSAMIterator pi = new ParallelSAMIterator(end1Reader, end2Reader);

        List<SAMRecord> end1 ;
        List<SAMRecord> end2 ;

        SAMRecord r1 = null, r2 = null; // alignments selected for ends 1 and 2, respectively, out of (possibly) multiple alternative placements

        while ( pi.hasNext() ) {

            fragments_seen++;

            Pair< List<SAMRecord>, List<SAMRecord> > ends = pi.next();

            end1 = ends.getFirst();
            end2 = ends.getSecond();

            if ( end1.size() == 0 ) {
                // nothing at all for end1, choose best from end2 and save
                end1_missing++;
                r1 = null;
                r2 = selectBestSingleEnd(end2);
                if ( AlignmentUtils.isReadUnmapped(r2) ) both_unmapped++;
                else end1_unmapped++;
 //               setPairingInformation(r1,r2);
 //               outWriter.addAlignment(r2);
                continue;
            }
            if ( end2.size() == 0 ) {
                // nothing at all for end2, choose best from end1 and save
                end2_missing++;
                r1 = selectBestSingleEnd(end1);
                r2 = null;
                if ( AlignmentUtils.isReadUnmapped(r1) ) both_unmapped++;
                else end2_unmapped++;
//                setPairingInformation(r1,r2);
//                outWriter.addAlignment(r1);
                continue;
            }

            if ( end1.size() == 1 && end2.size() == 1 ) {
                // unique alignments on both ends: not much to do, just save as a pair
                r1 = end1.get(0);
                r2 = end2.get(0);
                if ( AlignmentUtils.isReadUnmapped(r1) ) {
                    if ( AlignmentUtils.isReadUnmapped(r2) ) both_unmapped++;
                    else end1_unmapped++;
                } else {
                    if ( AlignmentUtils.isReadUnmapped(r2)) end2_unmapped++;
                    else {
                        // both mapped
                        both_mapped++;

                        if ( r1.getMappingQuality() >= MINQ.intValue() &&
                             r2.getMappingQuality() >= MINQ.intValue() ) {
                              both_unique++;
                              if ( r1.getReferenceIndex() == r2.getReferenceIndex() &&
                                   orientation(r1,r2) == PairOrientation.INNER ) {
                                  proper_pair++;
                                  setPairingInformation(r1,r2);
                                  outWriter.addAlignment(r1);
                                  outWriter.addAlignment(r2);
                              }
                        }
                    }
                }
//                setPairingInformation(r1,r2);
//                outWriter.addAlignment(r1);
//                outWriter.addAlignment(r2);
                continue;
            }

            if ( end1.size() == 1 && AlignmentUtils.isReadUnmapped(end1.get(0)) ) {
                // special case: multiple alignments for end2 but end1 is unmapped: just select best for end2
                r1 = end1.get(0);
                r2 = selectBestSingleEnd(end2);
                if ( AlignmentUtils.isReadUnmapped(r2) ) both_unmapped++;
                else end1_unmapped++;
//                setPairingInformation(r1,r2);
//                outWriter.addAlignment(r1);
//                outWriter.addAlignment(r2);
                continue;
            }

            if ( end2.size() == 1 && AlignmentUtils.isReadUnmapped(end2.get(0)) ) {
                // special case: multiple alignments for end1 but end2 is unmapped: just select best for end1
                r1 = selectBestSingleEnd(end1);
                r2 = end2.get(0);
                if ( AlignmentUtils.isReadUnmapped(r1) ) both_unmapped++;
                else end2_unmapped++;
 //               setPairingInformation(r1,r2);
 //               outWriter.addAlignment(r1);
 //               outWriter.addAlignment(r2);
                continue;
            }

            // ok, if we are here then we got both ends and multiple alignments in at least one end.
            // Let's loop through candidates and choose the best pair:
/*
            List<Pairing> good = new ArrayList<Pairing>();
            List<Pairing> bad = new ArrayList<Pairing>();
            double best_good = INFINITY;

            for ( SAMRecord candidate1 : end1 ) {
                for ( SAMRecord candidate2 : end2 ) {
                    if ( candidate1.getReferenceIndex() == candidate2.getReferenceIndex()
                         && orientation(candidate1,candidate2)==PairOrientation.INNER ) {
                        double score = pairingScore(candidate1, candidate2);
                    }
                }
            }
*/
            r1 = selectUniqueSingleEnd(end1,MINQ.intValue());
            r2 = selectUniqueSingleEnd(end2,MINQ.intValue());
            if ( r1 == null || r2 == null ) continue;

            both_mapped++;
            both_unique++;
            if ( r1.getReferenceIndex() == r2.getReferenceIndex() &&
                 orientation(r1,r2) == PairOrientation.INNER ) {
                      proper_pair++;
                      setPairingInformation(r1,r2);
                      outWriter.addAlignment(r1);
                      outWriter.addAlignment(r2);
            }


        }

        pi.close();
        outWriter.close();

        System.out.println("Pairs with end1 missing: "+end1_missing);
        System.out.println("Pairs with end2 missing: "+end2_missing);
        System.out.println("Pairs with only end1 unmapped (including missing): "+end1_unmapped);
        System.out.println("Pairs with only end2 unmapped (including missing): "+end2_unmapped);
        System.out.println("Pairs with both ends unmapped (including one missing): "+both_unmapped);
        System.out.println("Pairs with both ends mapped: "+both_mapped);
        System.out.println("Pairs with both ends mapped uniquely (based on MINQ="+MINQ+"): "+both_unique);
        System.out.println("Pairs with both ends mapped uniquely and properly (written into output): "+proper_pair);
        return 0;
    }

    private Query<Pair<SAMRecord,SAMRecord> > parseConditions(String filter) {

        filter = filter.trim();
        Query<Pair<SAMRecord,SAMRecord>> result1, result2;

        int level = 0; // parentheses level

        for ( int i = 0 ; i < filter.length() ; i++ ) {
            switch ( filter.charAt(i) ) {
            case '(': level++; break;
            case ')': level--;
                      if ( level < 0 ) throw new RuntimeException("Too many closing parentheses in the expression.");
                      break;
            case '&': if ( level > 0 ) break; // parenthised expression - not now!
                        // we are at level 0: parse expressions to the left and to the right of '&' operator
                      return new CompositeQuery<Pair<SAMRecord,SAMRecord> >(parseConditions(filter.substring(0,i)),
                                                                        parseConditions(filter.substring(i+1)),
                                                                        Query.Operator.AND);
            case '|': if ( level > 0 ) break; // inside parenthised expression - keep scanning, we'll process the whole expression later
                          // we are at level 0: parse expressions to the left and to the right of '|' operator
                      return new CompositeQuery<Pair<SAMRecord,SAMRecord> >(parseConditions(filter.substring(0,i)),
                                                                            parseConditions(filter.substring(i+1)),
                                                                            Query.Operator.OR);
            default: break;
            }
        }

        if ( level > 0 ) throw new RuntimeException("Too many opening parentheses in the expression.");

        // if we ended up here, this is either a single parenthized expression or a primitive.
        // filter was trimmed; if it is a parenthized expression, ( and ) should be first/last symbols:
        if ( filter.charAt(0) == '(' && filter.charAt(filter.length()-1) == ')')
            return parseConditions(filter.substring(1,filter.length()-1));

        // ok, it's a primitive:
        int equal_pos = filter.indexOf('=');
        if ( equal_pos < 0  ) { // it's not a <tag>=<value> expression, but a logical primitive
            if ( "PAIRSONLY".equals(filter) ) return new BothEndsMappedQuery();
        }

        return null;
    }

    /**
     * Utility method: checks if the two headers are the same. Returns the first one if they are,
     * a non-NULL one if the other one is NULL, or NULL if both headers are NULL. If headers are
     * both not NULL and are not the same, a RuntimeException is thrown. 
     * @param h1
     * @param h2
     * @return true if the headers are the same
     */
    private SAMFileHeader checkHeaders(SAMFileHeader h1, SAMFileHeader h2) {

        if ( h1 == null ) return h2;
        if ( h2 == null ) return h1;

//        if ( ! h1.getReadGroups().equals(h2.getReadGroups())) throw new RuntimeException("Read groups in the two input files do not match");
//        if ( ! h1.getSequenceDictionary().equals(h2.getSequenceDictionary()) ) throw new RuntimeException("Sequence dictionaries in the two input files do not match");
//        if ( ! h1.getProgramRecords().equals(h2.getProgramRecords()) ) throw new RuntimeException("Program records in the two input files do not match");
        if ( ! h1.equals(h2) ) throw new RuntimeException("Headers in the two input files do not match");
        return h1;
    }

    /** Given a list of alignments, returns the one with the best mapping quality score.
     *  If there is more than one alignment with the same score (got to be 0), one of
     *  these best-scoring alignments will be returned at random.
     * @param l
     * @return
     */
    private SAMRecord selectBestSingleEnd(List<SAMRecord> l) {
        if ( l.size() == 0 ) return null; // should not happen; just don not want to crash here, but somewhere else
        if ( l.size() == 1 ) return l.get(0); // not much choice...
        int best_qual = -1;
        int n_unmapped = 0;
        List<SAMRecord> best = new ArrayList<SAMRecord>();

        for ( SAMRecord r : l ) {
            if ( r.getReadUnmappedFlag() ) {
                // paranoid; if there are ANY alignments available, there should not be any "unmapped" records among them;
                // and if the read is "unmapped" indeed, then there should be no other alignments reported
                n_unmapped++;
                continue;
            }
            if ( r.getMappingQuality() > best_qual) {
                best_qual = r.getMappingQuality();
                best.clear();
                best.add(r);
                continue;
            }
            if ( r.getMappingQuality() == best_qual ) best.add(r);
        }
        if ( best.size() == 0 ) throw new RuntimeException("Currently Unsupported: SAM file might be not fully compliant. "+
                                                                        "Multiple 'unmapped' records found for read "+l.get(0).getReadName());
        if ( best.size() == 1 ) return best.get(0);
        if ( best_qual != 0 ) throw new RuntimeException("Multiple alignments for the same read found with non-zero score. "+
                                            "Read: "+l.get(0).getReadName()+" best score: "+best_qual);
        return best.get((int)(Math.random()*best.size()));
    }

    private SAMRecord selectUniqueSingleEnd(List<SAMRecord> l, int minq) {
        if ( l.size() == 0 ) return null; // should not happen; just don not want to crash here, but somewhere else
        if ( l.size() == 1 ) {
            if ( l.get(0).getMappingQuality() >= minq ) return l.get(0);
            else return null; // not unique enough
        }

        int n_unmapped = 0;
        List<SAMRecord> best = new ArrayList<SAMRecord>();

        for ( SAMRecord r : l ) {
            if ( AlignmentUtils.isReadUnmapped(r) ) {
                // paranoid; if there are ANY alignments available, there should not be any "unmapped" records among them;
                // and if the read is "unmapped" indeed, then there should be no other alignments reported
                n_unmapped++;
                continue;
            }
            if ( r.getMappingQuality() >= minq ) {
                best.add(r);
                continue;
            }
        }
        if ( best.size() == 0 ) return null; // no unique alignment
        if ( best.size() > 1 ) {
            for ( SAMRecord r : best ) {
                System.out.println("READ "+r.getReadName()+" mapQ="+r.getMappingQuality()+" at="+r.getReferenceName()+
                        ":"+r.getAlignmentStart()+"("+(r.getReadNegativeStrandFlag()?"-":"+")+") cig="+r.getCigarString());
            }
            throw new RuntimeException("Multiple alignments for read "+l.get(0).getReadName()+", all with Q>="+minq);
        }

        return best.get(0);
    }

    /**
     * Assumes that alignments r1 and r2 are the two ends of a selected mate pair, and sets pairing flags and reciprocal mate
     * mapping values for each of them (so that e.g. r2.getMateAlignmentStart() is properly set to
     * r1.getAlignmentStart(), etc). Any one of the two alignments can be null, in which case it is assumed
     * that the corresponding end is unmapped, and the flags/values in the other end will be set accordingly.
     * If both r1 and r2 are null, an exception will be thrown.
     * @param r1 first end in a mate pair
     * @param r2 second end in a mate pair
     */
    private void setPairingInformation(SAMRecord r1, SAMRecord r2) {
        // set mate information (note that r1 and r2 can not be 'null' simultaneously):

        if ( r1 == null && r2 == null ) throw new RuntimeException("Both ends of the mate pair are passed as 'null'");

        if ( r1 != null ) {
            r1.setMateUnmappedFlag( r2==null ? true : AlignmentUtils.isReadUnmapped(r2) );
            r1.setMateReferenceIndex( r2==null ? SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX : r2.getReferenceIndex() );
            r1.setMateNegativeStrandFlag(r2==null ? false : r2.getReadNegativeStrandFlag() );
            r1.setMateAlignmentStart( r2 == null ? SAMRecord.NO_ALIGNMENT_START : r2.getAlignmentStart());

            r1.setFirstOfPairFlag(true);
            r1.setSecondOfPairFlag(false);
            r1.setReadPairedFlag(true);
        }
        if ( r2 != null ) {
            r2.setMateUnmappedFlag( r1==null ? true : AlignmentUtils.isReadUnmapped(r1) );
            r2.setMateReferenceIndex( r1==null ? SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX : r1.getReferenceIndex() );
            r2.setMateNegativeStrandFlag(r1==null ? false : r1.getReadNegativeStrandFlag() );
            r2.setMateAlignmentStart( r1 == null ? SAMRecord.NO_ALIGNMENT_START : r1.getAlignmentStart());
            r2.setFirstOfPairFlag(false);
            r2.setSecondOfPairFlag(true);
            r2.setReadPairedFlag(true);
        }
    }

    /**
     * Returns fragment length inferred from the two alignments, or 1,000,000,000 if reads align
     * on different chromosomes or are not mapped at all. NOTE: the returned fragment length is
     * start+read length of the rightmost alignment minus start of the leftmost one; it makes
     * sense only for "proper" pairs (leftmost forward, rightmost reverse); for all other pairs
     * the returned number does not allow for any meaningful interpretation. 
     * @param r1
     * @param r2
     * @return
     */
    private int fragmentSize(final SAMRecord r1, final SAMRecord r2) {
        if ( r1 == null || AlignmentUtils.isReadUnmapped(r1) ||
             r2 == null || AlignmentUtils.isReadUnmapped(r2) ||
                r1.getReferenceIndex() != r2.getReferenceIndex() ) return INFINITY;
        if ( r1.getAlignmentStart() <= r2.getAlignmentStart() )
             return ( r2.getAlignmentStart() + r2.getReadLength() - r1.getAlignmentStart());
        else return ( r1.getAlignmentStart() + r1.getReadLength() - r2.getAlignmentStart());
    }

    enum PairOrientation {
        INNER, OUTER, LEFT, RIGHT, NONE
    }

    /**
     * Returns orientation of the pair: INNER for "-->   <--", OUTER for "<--  -->",
     * LEFT for "<--   <--" and RIGHT for "-->    -->" (regardless of which read in a pair, 1 or 2,
     * actually maps to the left and which to the right). If any of the reads is null or unmapped, returns NONE.
     * If reads are on different contigs, they are still ordered according to the underlying contig order (by reference
     * index), and the returned value reflects their relative orientation as described above (however it does not seem
     * to be very useful in that case).
     * @param r1
     * @param r2
     * @return
     */
    private PairOrientation orientation(SAMRecord r1, SAMRecord r2) {

        if ( r1 == null || r2 == null || AlignmentUtils.isReadUnmapped(r1) || AlignmentUtils.isReadUnmapped(r2))
            return PairOrientation.NONE;

        SAMRecord left, right;

        if ( r1.getReferenceIndex() == r2.getReferenceIndex() ) {
            if ( r1.getAlignmentStart() <= r2.getAlignmentStart() ) {
                left = r1;
                right = r2;
            } else {
                left = r2;
                right = r1;
            }
        } else {
            if ( r1.getReferenceIndex() < r2.getReferenceIndex() ) {
                left = r1;
                right = r2;
            } else {
                left = r2;
                right = r1;
            }
        }

        if ( !  left.getReadNegativeStrandFlag() ) { // left is forward
            if ( right.getReadNegativeStrandFlag() ) return PairOrientation.INNER; // left is forward, right is reverse
            else return PairOrientation.RIGHT;  // left is forward, right is forward
        }  else {  // left is reverse
            if ( right.getReadNegativeStrandFlag() ) return PairOrientation.LEFT; // left is reverse, right is reverse
            else return PairOrientation.OUTER;  // left is reverse, right is forward
        }
    }

    class Pairing {
        SAMRecord r1;
        SAMRecord r2;
        double score;

        Pairing() {
            this(null,null,INFINITY);
        }

        Pairing(SAMRecord r1, SAMRecord r2) {
            this(r1,r2,INFINITY);
        }
        Pairing(SAMRecord r1, SAMRecord r2, double score) {
            this.r1 = r1;
            this.r2 = r2;
            this.score = score;
        }

        SAMRecord getFirst() { return r1; }
        SAMRecord getSecond() { return r2; }
        double getScore() { return score; }

        void setFirst(SAMRecord r) { r1 = r; }
        void setSecond(SAMRecord r) { r2 = r; }
        void setScore(double score) { this.score = score; }
    }

    private double pairingScore(final SAMRecord r1, final SAMRecord r2) {
        
        return ( ( r1.getMappingQuality() + r2.getMappingQuality() ) * Math.exp(1))  ;
    }

    interface Query<T> {
        boolean isSatisfied(T record) ;
        enum Operator { OR, AND };
    }

    class CompositeQuery<T> implements Query<T> {
        private Query<T> q1;
        private Query<T> q2;
        private Query.Operator type ; // 1 for 'and', 0 for 'or'

        CompositeQuery(Query<T> q1, Query<T> q2, Query.Operator type) {
            this.q1 = q1;
            this.q2 = q2;
            this.type = type;
        }

        public boolean isSatisfied(T record) {
            switch ( type ) {
            case AND: return q1.isSatisfied(record) && q2.isSatisfied(record);
            case OR:  return q1.isSatisfied(record) || q2.isSatisfied(record);
            default: throw new IllegalStateException("Unknown composite query operator");
            }
        }
    }

    class BothEndsMappedQuery implements Query< Pair<SAMRecord,SAMRecord> > {
        public boolean isSatisfied(Pair<SAMRecord,SAMRecord> p) {
            return ( ! AlignmentUtils.isReadUnmapped(p.getFirst()) && ! AlignmentUtils.isReadUnmapped(p.getSecond())) ;
        }
    }

}
