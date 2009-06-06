package org.broadinstitute.sting.playground.indels;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.picard.reference.ReferenceSequence;
import org.broadinstitute.sting.playground.utils.CountedObject;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 25, 2009
 * Time: 12:15:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentUtils {


    /** Returns number of mismatches in the alignment <code>r</code> to the reference sequence
     * <code>refSeq</code>. It is assumed that
     * the alignment starts at (1-based) position r.getAlignmentStart() on the specified, and all single-base mismatches
     * are counted in the alignment segments where both sequences are present. Insertions/deletions are skipped and do
     * not contribute to the error count returned by this method. 
     * @param r aligned read
     * @param refSeq reference sequence
     * @return number of single-base mismatches in the aligned segments (gaps on either of the sequences are skipped)
     */
    public static int numMismatches(SAMRecord r, ReferenceSequence refSeq) {
        byte[] ref = refSeq.getBases();
        if ( r.getReadUnmappedFlag() ) return 1000000;
        int i_ref = r.getAlignmentStart()-1; // position on the ref
        int i_read = 0;                    // position on the read
        int mm_count = 0; // number of mismatches
        Cigar c = r.getCigar();
        for ( int k = 0 ; k < c.numCigarElements() ; k++ ) {
            CigarElement ce = c.getCigarElement(k);
            switch( ce.getOperator() ) {
                case M:
                    for ( int l = 0 ; l < ce.getLength() ; l++, i_ref++, i_read++ ) {
                        if ( Character.toUpperCase(r.getReadString().charAt(i_read) ) == 'N' ) continue; // do not count N's ?
                        if ( Character.toUpperCase(r.getReadString().charAt(i_read) ) !=
                             Character.toUpperCase((char)ref[i_ref]) ) mm_count++;
                    }
                    break;
                case I: i_read += ce.getLength(); break;
                case D: i_ref += ce.getLength(); break;
                default: throw new RuntimeException("Unrecognized cigar element");
            }

        }
        return mm_count;
    }

    /**
     * mhanna - 11 May 2009 - stubbed out competing method that works with partial references. 
     * Computes number of mismatches in the read alignment to the refence <code>ref</code>
     * specified in the record <code>r</code>. Indels are completely <i>ignored</i> by this method:
     * only base mismatches in the alignment segments where both sequences are present are counted.
     * @param r
     * @return
     */
    public static int numMismatches(SAMRecord r, char[] ref) {
        if ( r.getReadUnmappedFlag() ) return 1000000;
        int i_ref = 0; // position on the ref
        int i_read = 0;                    // position on the read
        int mm_count = 0; // number of mismatches
        Cigar c = r.getCigar();
        for ( int k = 0 ; k < c.numCigarElements() ; k++ ) {
            CigarElement ce = c.getCigarElement(k);
            switch( ce.getOperator() ) {
                case M:
                    for ( int l = 0 ; l < ce.getLength() ; l++, i_ref++, i_read++ ) {
                        if ( Character.toUpperCase(r.getReadString().charAt(i_read) ) == 'N' ) continue; // do not count N's ?
                        if ( Character.toUpperCase(r.getReadString().charAt(i_read) ) !=
                             Character.toUpperCase(ref[i_ref]) ) mm_count++;
                    }
                    break;
                case I: i_read += ce.getLength(); break;
                case D: i_ref += ce.getLength(); break;
                default: throw new RuntimeException("Unrecognized cigar element");
            }

        }
        return mm_count;
    }

    // IMPORTANT NOTE: ALTHOUGH THIS METHOD IS EXTREMELY SIMILAR TO THE ONE ABOVE, WE NEED
    // TWO SEPARATE IMPLEMENTATIONS IN ORDER TO PREVENT JAVA STRINGS FROM FORCING US TO
    // PERFORM EXPENSIVE ARRAY COPYING WHEN TRYING TO GET A BYTE ARRAY...
    /** See {@link #numMismatches(SAMRecord, ReferenceSequence)}. This method implements same functionality
     * for reference sequence specified as conventional java string (of bases). By default, it is assumed that
     * the alignment starts at (1-based) position r.getAlignmentStart() on the reference <code>refSeq</code>.
     * See {@link #numMismatches(SAMRecord, String, int)} if this is not the case.
     */
    public static int numMismatches(SAMRecord r, String refSeq ) {
        if ( r.getReadUnmappedFlag() ) return 1000000;
        return numMismatches(r, refSeq, r.getAlignmentStart()-1);
     }

    /** Returns number of mismatches in the alignment <code>r</code> to the reference sequence
     * <code>refSeq</code> assuming the alignment starts at (ZERO-based) position <code>refIndex</code> on the
     * specified reference sequence; in other words, <code>refIndex</code> is used in place of alignment's own
     * getAlignmentStart() coordinate and the latter is never used. However, the structure of the alignment <code>r</code>
     * (i.e. it's cigar string with all the insertions/deletions it may specify) is fully respected.
     * 
     * @param r alignment
     * @param refSeq chunk of reference sequence that subsumes the alignment completely (if alignment runs out of 
     *                  the reference string, IndexOutOfBound exception will be thrown at runtime).
     * @param refIndex zero-based position, at which the alignment starts on the specified reference string. 
     * @return
     */
    public static int numMismatches(SAMRecord r, String refSeq, int refIndex) {
        int readIdx = 0;
        int mismatches = 0;
        String readSeq = r.getReadString();
        Cigar c = r.getCigar();
        for (int i = 0 ; i < c.numCigarElements() ; i++) {
            CigarElement ce = c.getCigarElement(i);
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0 ; j < ce.getLength() ; j++, refIndex++, readIdx++ ) {
                        if ( Character.toUpperCase(readSeq.charAt(readIdx)) != Character.toUpperCase(refSeq.charAt(refIndex)) )
                            mismatches++;
                    }
                    break;
                case I:
                    readIdx += ce.getLength();
                    break;
                case D:
                    refIndex += ce.getLength();
                    break;
                default: throw new StingException("Only M,I,D cigar elements are currently supported");
            }

        }
        return mismatches;
    }

    
    /** Reads through the alignment cigar and returns all indels found in the alignment as a collection
     * of Indel objects.
     * @param c
     * @param start
     * @return
     */
    public static Collection<Indel> extractIndels(final Cigar c, final int start) {
        //
        // firstpos,lastpos span of the indel will be interpreted as follows:
        // any alignment that ends strictly before firstpos or starts strictly after lastpos
        // on the *reference* (not inclusive!) does not overlap with an indel; in the case of
        // insertion it will result in firstpos > lastpos!
        //         lastpos
        //         |   firstpos
        //         |   |
        //         v   v
        // ---------III----- Ref  Insertion: bases I are not in the ref; any alignment that starts
        //                        after lastpos or ends before firstpos *on the reference*
        //                        is completely over the reference bases to the right or to
        //                        the left, respectively, of the insertion site
        //
        //      firstpos
        //      | lastpos
        //      | |
        //      v v
        //------------------ Ref   Deletion: any alignment that ends before firstpos or starts after lastpos
        // -----DDD--- alignment   on the reference does not overlap with the deletion
        int runninglength = start; // position on the original reference; start = alignment start position

        List<Indel> indels = new ArrayList<Indel>(4);

        if ( c.numCigarElements() <= 1 ) return indels; // most of the reads have no indels, save a few cycles by returning early

        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {

            final CigarElement ce = c.getCigarElement(i);
            Indel curr_indel = null;

            switch(ce.getOperator()) {
            case I:
                    curr_indel = new Indel(runninglength, ce.getLength(), Indel.IndelType.I);
                    if ( i == 0 ) System.out.println("WARNING: Indel at start!");
                    if ( i == c.numCigarElements() - 1) System.out.println("WARNING: Indel at end!");
                    break;
            case D: curr_indel = new Indel(runninglength, ce.getLength(), Indel.IndelType.D);
                    if ( i == 0 ) System.out.println("WARNING: Indel at start!");
                    if ( i == c.numCigarElements() - 1) System.out.println("WARNING: Indel at end!");
                    runninglength += ce.getLength();
                    break;
            case M: runninglength += ce.getLength(); break; // advance along the gapless block in the alignment
            default :
                throw new IllegalArgumentException("Unexpected operator in cigar string");
            }

            if ( curr_indel == null ) continue; // element was not an indel, go grab next element...

            indels.add(curr_indel); // this is a new indel. Add it.
        } // end for loop over all alignment cigar elements

        return indels;
    } // end extractIndels() method

    /** Reads through the alignment specified in the record and returns all indels found in the alignment as a collection
     * of Indel objects. If read is not mapped, silently returns an empty collection.
     * @param r
     * @return
     */
    public static Collection<Indel> extractIndels(SAMRecord r) {
        if ( r.getReadUnmappedFlag() ) return new ArrayList<Indel>();
        return extractIndels(r.getCigar(), r.getAlignmentStart());
    }

    /** Extracts indels from the specified record (@see #extractIndels()) and updates the provided tree object.
     * Multiple occurences of the same indel (i.e. with support from multiple reads) will be grouped together
     * in one counted object held by the set.
     *
     * @param r
     * @param t
     */
    public static void collectAndCountIndels(SAMRecord r, TreeSet<CountedObject<Indel> > t) {
        Collection<Indel> indels = AlignmentUtils.extractIndels(r);
        for ( Indel ind : indels ) {
            CountedObject<Indel> ci = new CountedObject<Indel>(ind);
            CountedObject<Indel> found = t.floor(ci);
//            CountedObject<Indel> found2 = t.ceiling(ci);

            if ( ci.equals( found ) ) found.increment(); // we did find our indel, advance the counter
            else t.add(ci); // this is a new indel. Add it.

        }

    }

    public static String toString(Cigar cig) {
        StringBuilder b = new StringBuilder();

        for ( int i = 0 ; i < cig.numCigarElements() ; i++ ) {
            char c='?';
            switch ( cig.getCigarElement(i).getOperator() ) {
            case M : c = 'M'; break;
            case D : c = 'D'; break;
            case I : c = 'I'; break;
            }
            b.append(c);
            b.append(cig.getCigarElement(i).getLength());
        }
        return b.toString();
    }
}
