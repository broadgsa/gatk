package org.broadinstitute.sting.playground.indels;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import edu.mit.broad.picard.reference.ReferenceSequence;
import org.broadinstitute.sting.playground.utils.CountedObject;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 25, 2009
 * Time: 12:15:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentUtils {

    /** Computes number of mismatches in the read alignment to the refence <code>ref</code>
     * specified in the record <code>r</code>. Indels are completely <i>ignored</i> by this method:
     * only base mismatches in the alignment segments where both sequences are present are counted.
     * @param r
     * @param ref
     * @return
     */
    public static int numMismatches(SAMRecord r, ReferenceSequence ref) {
        return numMismatches(r, ref.getBases().toString());
    }

    public static int numMismatches(SAMRecord r, String ref) {
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
                             Character.toUpperCase(ref.charAt(i_ref)) ) mm_count++;
                    }
                    break;
                case I: i_read += ce.getLength(); break;
                case D: i_ref += ce.getLength(); break;
                default: throw new RuntimeException("Unrecognized cigar element");
            }

        }
        return mm_count;
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
