package org.broadinstitute.sting.utils;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.pileup.*;


public class AlignmentUtils {

    private static class MismatchCount {
        public int numMismatches = 0;
        public long mismatchQualities = 0;
    }

    /** Returns number of mismatches in the alignment <code>r</code> to the reference sequence
     * <code>refSeq</code> assuming the alignment starts at (ZERO-based) position <code>refIndex</code> on the
     * specified reference sequence; in other words, <code>refIndex</code> is used in place of alignment's own
     * getAlignmentStart() coordinate and the latter is never used. However, the structure of the alignment <code>r</code>
     * (i.e. it's cigar string with all the insertions/deletions it may specify) is fully respected.
     *
     * THIS CODE ASSUMES THAT ALL BYTES COME FROM UPPERCASED CHARS.
     * 
     * @param r alignment
     * @param refSeq chunk of reference sequence that subsumes the alignment completely (if alignment runs out of 
     *                  the reference string, IndexOutOfBound exception will be thrown at runtime).
     * @param refIndex zero-based position, at which the alignment starts on the specified reference string. 
     * @return the number of mismatches
     */
    public static int numMismatches(SAMRecord r, byte[] refSeq, int refIndex) {
        return getMismatchCount(r, refSeq, refIndex).numMismatches;
    }

    public static int numMismatches(SAMRecord r, String refSeq, int refIndex ) {
        if ( r.getReadUnmappedFlag() ) return 1000000;
        return numMismatches(r, StringUtil.stringToBytes(refSeq), refIndex);
     }

    public static long mismatchingQualities(SAMRecord r, byte[] refSeq, int refIndex) {
        return getMismatchCount(r, refSeq, refIndex).mismatchQualities;
    }

    public static long mismatchingQualities(SAMRecord r, String refSeq, int refIndex ) {
        if ( r.getReadUnmappedFlag() ) return 1000000;
        return numMismatches(r, StringUtil.stringToBytes(refSeq), refIndex);
     }

    private static MismatchCount getMismatchCount(SAMRecord r, byte[] refSeq, int refIndex) {
        MismatchCount mc = new MismatchCount();

        int readIdx = 0;
        byte[] readSeq = r.getReadBases();
        Cigar c = r.getCigar();
        for (int i = 0 ; i < c.numCigarElements() ; i++) {
            CigarElement ce = c.getCigarElement(i);
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0 ; j < ce.getLength() ; j++, refIndex++, readIdx++ ) {
                        if ( refIndex >= refSeq.length )
                            continue;
                        byte refChr = refSeq[refIndex];
                        byte readChr = readSeq[readIdx];
                        // Note: we need to count X/N's as mismatches because that's what SAM requires
                        //if ( BaseUtils.simpleBaseToBaseIndex(readChr) == -1 ||
                        //     BaseUtils.simpleBaseToBaseIndex(refChr)  == -1 )
                        //    continue; // do not count Ns/Xs/etc ?
                        if ( readChr != refChr ) {
                            mc.numMismatches++;
                            mc.mismatchQualities += r.getBaseQualities()[readIdx];
                        }
                    }
                    break;
                case I:
                case S:
                    readIdx += ce.getLength();
                    break;
                case D:
                case N:
                    refIndex += ce.getLength();
                    break;
                default: throw new StingException("The " + ce.getOperator() + " cigar element is not currently supported");
            }

        }
        return mc;
    }

    /** Returns the number of mismatches in the pileup within the given reference context.
     *
     * @param pileup  the pileup with reads
     * @param ref     the reference context
     * @param ignoreTargetSite     if true, ignore mismatches at the target locus (i.e. the center of the window)
     * @return the number of mismatches
     */
    public static int mismatchesInRefWindow(ReadBackedPileup pileup, ReferenceContext ref, boolean ignoreTargetSite) {
        int mismatches = 0;
        for ( PileupElement p : pileup )
            mismatches += mismatchesInRefWindow(p, ref, ignoreTargetSite);
        return mismatches;
    }

    /** Returns the number of mismatches in the pileup element within the given reference context.
     *
     * @param p       the pileup element
     * @param ref     the reference context
     * @param ignoreTargetSite     if true, ignore mismatches at the target locus (i.e. the center of the window)
     * @return the number of mismatches
     */
    public static int mismatchesInRefWindow(PileupElement p, ReferenceContext ref, boolean ignoreTargetSite) {
        return mismatchesInRefWindow(p, ref, ignoreTargetSite, false);
    }

    /** Returns the number of mismatches in the pileup element within the given reference context.
     *
     * @param p       the pileup element
     * @param ref     the reference context
     * @param ignoreTargetSite     if true, ignore mismatches at the target locus (i.e. the center of the window)
     * @param qualitySumInsteadOfMismatchCount if true, return the quality score sum of the mismatches rather than the count
     * @return the number of mismatches
     */
    public static int mismatchesInRefWindow(PileupElement p, ReferenceContext ref, boolean ignoreTargetSite, boolean qualitySumInsteadOfMismatchCount) {
        int sum = 0;

        int windowStart = (int)ref.getWindow().getStart();
        int windowStop = (int)ref.getWindow().getStop();
        char[] refBases = ref.getBases();
        byte[] readBases = p.getRead().getReadBases();
        byte[] readQualities = p.getRead().getBaseQualities();
        Cigar c = p.getRead().getCigar();

        int readIndex = 0;
        int currentPos = p.getRead().getAlignmentStart();
        int refIndex = Math.max(0, currentPos - windowStart);

        for (int i = 0 ; i < c.numCigarElements() ; i++) {
            CigarElement ce = c.getCigarElement(i);
            int cigarElementLength = ce.getLength();
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0; j < cigarElementLength; j++, readIndex++, currentPos++) {
                        // are we past the ref window?
                        if ( currentPos > windowStop )
                            break;

                        // are we before the ref window?
                        if ( currentPos < windowStart )
                            continue;

                        char refChr = refBases[refIndex++];

                        // do we need to skip the target site?
                        if ( ignoreTargetSite && ref.getLocus().getStart() == currentPos )
                            continue;

                        char readChr = (char)readBases[readIndex];
                        if ( Character.toUpperCase(readChr) != Character.toUpperCase(refChr) )                       
                            sum += (qualitySumInsteadOfMismatchCount) ? readQualities[readIndex] : 1;
                    }
                    break;
                case I:
                case S:
                    readIndex += cigarElementLength;
                    break;
                case D:
                case N:
                    currentPos += cigarElementLength;
                    if ( currentPos > windowStart )
                        refIndex += Math.min(cigarElementLength, currentPos - windowStart);
                    break;
                default:
                    // fail silently
                    return 0;
            }
        }

        return sum;
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

    public static char[] alignmentToCharArray( final Cigar cigar, final char[] read, final char[] ref ) {

        final char[] alignment = new char[read.length];
        int refPos = 0;
        int alignPos = 0;

        for ( int iii = 0 ; iii < cigar.numCigarElements() ; iii++ ) {

            final CigarElement ce = cigar.getCigarElement(iii);

            switch( ce.getOperator() ) {
            case I:
            case S:
                for ( int jjj = 0 ; jjj < ce.getLength(); jjj++ ) {
                    alignment[alignPos++] = '+';
                }
                break;
            case D:
            case N:
                refPos++;
                break;
            case M:
                for ( int jjj = 0 ; jjj < ce.getLength(); jjj++ ) {
                    alignment[alignPos] = ref[refPos];
                    alignPos++;
                    refPos++;
                }
                break;
            default:
                throw new StingException( "Unsupported cigar operator: " + ce.getOperator() );
            }
        }
        return alignment;
    }

    /**
     * Due to (unfortunate) multiple ways to indicate that read is unmapped allowed by SAM format
     * specification, one may need this convenience shortcut. Checks both 'read unmapped' flag and
     * alignment reference index/start.
     * @param r record
     * @return true if read is unmapped
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
