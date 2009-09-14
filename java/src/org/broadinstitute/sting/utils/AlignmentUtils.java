package org.broadinstitute.sting.utils;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.picard.reference.ReferenceSequence;
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
                        char refChr = (char)ref[i_ref];
                        char readChr = r.getReadString().charAt(i_read);
                        if ( BaseUtils.simpleBaseToBaseIndex(readChr) == -1 ||
                             BaseUtils.simpleBaseToBaseIndex(refChr)  == -1 )
                            continue; // do not count Ns/Xs/etc ?
                        if ( Character.toUpperCase(readChr) != Character.toUpperCase(refChr) )
                            mm_count++;
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
                        char refChr = (char)ref[i_ref];
                        char readChr = r.getReadString().charAt(i_read);
                        if ( BaseUtils.simpleBaseToBaseIndex(readChr) == -1 ||
                             BaseUtils.simpleBaseToBaseIndex(refChr)  == -1 )
                            continue; // do not count Ns/Xs/etc ?
                        if ( Character.toUpperCase(readChr) != Character.toUpperCase(refChr) )
                             mm_count++;
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
                        if ( refIndex >= refSeq.length() )
                            continue;
                        char refChr = refSeq.charAt(refIndex);
                        char readChr = readSeq.charAt(readIdx);
                        // Note: we need to count X/N's as mismatches because that's what SAM requires
                        //if ( BaseUtils.simpleBaseToBaseIndex(readChr) == -1 ||
                        //     BaseUtils.simpleBaseToBaseIndex(refChr)  == -1 )
                        //    continue; // do not count Ns/Xs/etc ?
                        if ( Character.toUpperCase(readChr) != Character.toUpperCase(refChr) )
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

    /** Returns number of alignment blocks (continuous stretches of aligned bases) in the specified alignment.
     * This method follows closely the SAMRecord::getAlignmentBlocks() implemented in samtools library, but
     * it only counts blocks without actually allocating and filling the list of blocks themselves. Hence, this method is
     * a much more efficient alternative to r.getAlignmentBlocks.size() in the situations when this number is all that is needed.
     * Formally, this method simply returns the number of M elements in the cigar. 
     * @param r alignment
     * @return number of continuous alignment blocks (i.e. 'M' elements of the cigar; all indel and clipping elements are ignored).
     */
    public static int getNumAlignmentBlocks(final SAMRecord r) {
    	int n = 0;
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;
 
        for (final CigarElement e : cigar.getCigarElements()) {
        	if (e.getOperator() == CigarOperator.M ) n++;  
        }

    	return n;
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
            b.append(cig.getCigarElement(i).getLength());
            b.append(c);
        }
        return b.toString();
    }


    public static String alignmentToString(final Cigar cigar,final  String seq, final String ref, final int posOnRef ) {
        return alignmentToString( cigar, seq, ref, posOnRef, 0 );
    }

    public static String cigarToString(Cigar cig) {
        if ( cig == null )
            return "null";

        StringBuilder b = new StringBuilder();

        for ( int i = 0 ; i < cig.numCigarElements() ; i++ ) {
            char c='?';
            switch ( cig.getCigarElement(i).getOperator() ) {
                case M : c = 'M'; break;
                case D : c = 'D'; break;
                case I : c = 'I'; break;
            }
            b.append(cig.getCigarElement(i).getLength());
            b.append(c);
        }
        return b.toString();
    }

    public static String alignmentToString(final Cigar cigar,final  String seq, final String ref, final int posOnRef, final int posOnRead ) {
        int readPos = posOnRead;
        int refPos = posOnRef;
        
        StringBuilder refLine = new StringBuilder();
        StringBuilder readLine = new StringBuilder();

        for ( int i = 0 ; i < posOnRead ; i++ ) {
            refLine.append( ref.charAt( refPos - readPos + i ) );
            readLine.append( seq.charAt(i) ) ;
        }

        for ( int i = 0 ; i < cigar.numCigarElements() ; i++ ) {

            final CigarElement ce = cigar.getCigarElement(i);

            switch(ce.getOperator()) {
            case I:
                for ( int j = 0 ; j < ce.getLength(); j++ ) {
                    refLine.append('+');
                    readLine.append( seq.charAt( readPos++ ) );
                }
                break;
            case D:
                for ( int j = 0 ; j < ce.getLength(); j++ ) {
                    readLine.append('*');
                    refLine.append( ref.charAt( refPos++ ) );
                }
                break;
            case M:
                for ( int j = 0 ; j < ce.getLength(); j++ ) {
                    refLine.append(ref.charAt( refPos++ ) );
                    readLine.append( seq.charAt( readPos++ ) );
                }
                break;
            default: throw new StingException("Unsupported cigar operator: "+ce.getOperator() );
            }
        }
        refLine.append('\n');
        refLine.append(readLine);
        refLine.append('\n');
        return refLine.toString();
    }

    /**
     * Due to (unfortunate) multiple ways to indicate that read is unmapped allowed by SAM format
     * specification, one may need this convenience shortcut. Checks both 'read unmapped' flag and
     * alignment reference index/start.
     * @param r
     * @return
     */
    public static boolean isReadUnmapped(final SAMRecord r) {
        if ( r.getReadUnmappedFlag() ) return true;

        // our life would be so much easier if all sam files followed the specs. In reality,
        // sam files (including those generated by maq or bwa) miss headers alltogether. When
        // reading such a SAM file, reference name is set, but since there is no sequence dictionary,
        // null is always returned for referenceIndex. Let's be paranoid here, and make sure that
        // we do not call the read "unmapped" when it has only reference name set with ref. index missing
        // or vice versa.
        if ( ( r.getReferenceIndex() != null && r.getReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX
                || r.getReferenceName() != null && r.getReferenceName() != SAMRecord.NO_ALIGNMENT_REFERENCE_NAME )
          &&  r.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START ) return false  ;
        return true;
    }
}
