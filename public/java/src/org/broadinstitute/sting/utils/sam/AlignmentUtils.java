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

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;


public class AlignmentUtils {

    public static class MismatchCount {
        public int numMismatches = 0;
        public long mismatchQualities = 0;
    }

    public static long mismatchingQualities(SAMRecord r, byte[] refSeq, int refIndex) {
        return getMismatchCount(r, refSeq, refIndex).mismatchQualities;
    }

    public static MismatchCount getMismatchCount(SAMRecord r, byte[] refSeq, int refIndex) {
        return getMismatchCount(r,refSeq,refIndex,0,r.getReadLength());
    }

    // todo -- this code and mismatchesInRefWindow should be combined and optimized into a single
    // todo -- high performance implementation.  We can do a lot better than this right now
    public static MismatchCount getMismatchCount(SAMRecord r, byte[] refSeq, int refIndex, int startOnRead, int nReadBases) {
        MismatchCount mc = new MismatchCount();

        int readIdx = 0;
        int endOnRead = startOnRead + nReadBases - 1; // index of the last base on read we want to count
        byte[] readSeq = r.getReadBases();
        Cigar c = r.getCigar();
        for (int i = 0 ; i < c.numCigarElements() ; i++) {

            if ( readIdx > endOnRead ) break;

            CigarElement ce = c.getCigarElement(i);
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0 ; j < ce.getLength() ; j++, refIndex++, readIdx++ ) {
                        if ( refIndex >= refSeq.length )
                            continue;
                        if ( readIdx < startOnRead ) continue;
                        if ( readIdx > endOnRead ) break;
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
                case H:
                case P:
                    break;
                default: throw new ReviewedStingException("The " + ce.getOperator() + " cigar element is not currently supported");
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

        int windowStart = ref.getWindow().getStart();
        int windowStop = ref.getWindow().getStop();
        byte[] refBases = ref.getBases();
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

                        byte refChr = refBases[refIndex++];

                        // do we need to skip the target site?
                        if ( ignoreTargetSite && ref.getLocus().getStart() == currentPos )
                            continue;

                        byte readChr = readBases[readIndex];
                        if ( readChr != refChr )
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
                case H:
                case P:
                    break;
            }
        }

        return sum;
    }

    /** Returns the number of mismatches in the pileup element within the given reference context.
     *
     * @param read          the SAMRecord
     * @param ref           the reference context
     * @param maxMismatches the maximum number of surrounding mismatches we tolerate to consider a base good
     * @param windowSize    window size (on each side) to test
     * @return a bitset representing which bases are good
     */
    public static BitSet mismatchesInRefWindow(SAMRecord read, ReferenceContext ref, int maxMismatches, int windowSize) {
        // first determine the positions with mismatches
        int readLength = read.getReadLength();
        BitSet mismatches = new BitSet(readLength);

        // it's possible we aren't starting at the beginning of a read,
        //  and we don't need to look at any of the previous context outside our window
        //  (although we do need future context)
        int readStartPos = Math.max(read.getAlignmentStart(), ref.getLocus().getStart() - windowSize);
        int currentReadPos = read.getAlignmentStart();

        byte[] refBases = ref.getBases();
        int refIndex = readStartPos - ref.getWindow().getStart();
        if ( refIndex < 0 ) {
            throw new IllegalStateException("When calculating mismatches, we somehow don't have enough previous reference context for read " + read.getReadName() + " at position " + ref.getLocus());
        }

        byte[] readBases = read.getReadBases();
        int readIndex = 0;

        Cigar c = read.getCigar();

        for (int i = 0 ; i < c.numCigarElements() ; i++) {
            CigarElement ce = c.getCigarElement(i);
            int cigarElementLength = ce.getLength();
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0; j < cigarElementLength; j++, readIndex++) {
                        // skip over unwanted bases
                        if ( currentReadPos++ < readStartPos )
                            continue;

                        // this is possible if reads extend beyond the contig end
                        if ( refIndex >= refBases.length )
                            break;

                        byte refChr = refBases[refIndex];
                        byte readChr = readBases[readIndex];
                        if ( readChr != refChr )
                            mismatches.set(readIndex);

                        refIndex++;
                    }
                    break;
                case I:
                case S:
                    readIndex += cigarElementLength;
                    break;
                case D:
                case N:
                    if ( currentReadPos >= readStartPos )
                        refIndex += cigarElementLength;
                    currentReadPos += cigarElementLength;
                    break;
                case H:
                case P:
                    break;
            }
        }

        // all bits are set to false by default
        BitSet result = new BitSet(readLength);

        int currentPos = 0, leftPos = 0, rightPos;
        int mismatchCount = 0;

        // calculate how many mismatches exist in the windows to the left/right
        for ( rightPos = 1; rightPos <= windowSize && rightPos < readLength; rightPos++) {
            if ( mismatches.get(rightPos) )
                mismatchCount++;
        }
        if ( mismatchCount <= maxMismatches )
            result.set(currentPos);

        // now, traverse over the read positions
        while ( currentPos < readLength ) {
            // add a new rightmost position
            if ( rightPos < readLength && mismatches.get(rightPos++) )
                mismatchCount++;
            // re-penalize the previous position
            if ( mismatches.get(currentPos++) )
                mismatchCount++;
            // don't penalize the current position
            if ( mismatches.get(currentPos) )
                mismatchCount--;
            // subtract the leftmost position
            if ( leftPos < currentPos - windowSize && mismatches.get(leftPos++) )
                    mismatchCount--;

            if ( mismatchCount <= maxMismatches )
                result.set(currentPos);
        }

        return result;
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

    public static int getNumAlignedBases(final SAMRecord r) {
    	int n = 0;
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;

        for (final CigarElement e : cigar.getCigarElements()) {
        	if (e.getOperator() == CigarOperator.M ) { n += e.getLength(); }
        }

    	return n;
    }

    public static byte[] alignmentToByteArray( final Cigar cigar, final byte[] read, final byte[] ref ) {

        final byte[] alignment = new byte[read.length];
        int refPos = 0;
        int alignPos = 0;

        for ( int iii = 0 ; iii < cigar.numCigarElements() ; iii++ ) {

            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch( ce.getOperator() ) {
            case I:
            case S:
                for ( int jjj = 0; jjj < elementLength; jjj++ ) {
                    alignment[alignPos++] = '+';
                }
                break;
            case D:
            case N:
                refPos += elementLength;
                break;
            case M:
                for ( int jjj = 0; jjj < elementLength; jjj++ ) {
                    alignment[alignPos] = ref[refPos];
                    alignPos++;
                    refPos++;
                }
                break;
            case H:
            case P:
                break;
            default:
                throw new ReviewedStingException( "Unsupported cigar operator: " + ce.getOperator() );
            }
        }
        return alignment;
    }

    public static int calcAlignmentByteArrayOffset( final Cigar cigar, int pileupOffset, final int alignmentStart, final int refLocus ) {

        boolean atDeletion = false;
        if(pileupOffset == -1) {
            atDeletion = true;
            pileupOffset = refLocus - alignmentStart;
            final CigarElement ce = cigar.getCigarElement(0);
            if( ce.getOperator() == CigarOperator.S ) {
                pileupOffset += ce.getLength();
            }
        }
       int pos = 0;
       int alignmentPos = 0;
       for ( int iii = 0 ; iii < cigar.numCigarElements() ; iii++ ) {

            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch( ce.getOperator() ) {
            case I:
            case S:
                pos += elementLength;
                if( pos >= pileupOffset ) {
                    return alignmentPos;
                }
                break;
            case D:
            case N:
                if(!atDeletion) {
                    alignmentPos += elementLength;
                } else {
                    if( pos + elementLength - 1 >= pileupOffset ) {
                        return alignmentPos + (pileupOffset - pos);
                    } else {
                        pos += elementLength;
                        alignmentPos += elementLength;
                    }
                }
                break;
            case M:
                if( pos + elementLength - 1 >= pileupOffset ) {
                    return alignmentPos + (pileupOffset - pos);
                } else {
                    pos += elementLength;
                    alignmentPos += elementLength;
                }
                break;
            case H:
            case P:
                break;
            default:
                throw new ReviewedStingException( "Unsupported cigar operator: " + ce.getOperator() );
            }
        }
        return alignmentPos;
    }

   public static byte[] readToAlignmentByteArray( final Cigar cigar, final byte[] read ) {

       int alignmentLength = 0;
       for ( int iii = 0 ; iii < cigar.numCigarElements() ; iii++ ) {

            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch( ce.getOperator() ) {
            case I:
            case S:
                break;
            case D:
            case N:
                alignmentLength += elementLength;
                break;
            case M:
                alignmentLength += elementLength;
                break;
            case H:
            case P:
                break;
            default:
                throw new ReviewedStingException( "Unsupported cigar operator: " + ce.getOperator() );
            }
        }

        final byte[] alignment = new byte[alignmentLength];
        int alignPos = 0;
        int readPos = 0;
        for ( int iii = 0 ; iii < cigar.numCigarElements() ; iii++ ) {

            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch( ce.getOperator() ) {
            case I:
                if( alignPos > 0 ) {
                    if( alignment[alignPos-1] == BaseUtils.A ) { alignment[alignPos-1] = PileupElement.A_FOLLOWED_BY_INSERTION_BASE; }
                    else if( alignment[alignPos-1] == BaseUtils.C ) { alignment[alignPos-1] = PileupElement.C_FOLLOWED_BY_INSERTION_BASE; }
                    else if( alignment[alignPos-1] == BaseUtils.T ) { alignment[alignPos-1] = PileupElement.T_FOLLOWED_BY_INSERTION_BASE; }
                    else if( alignment[alignPos-1] == BaseUtils.G ) { alignment[alignPos-1] = PileupElement.G_FOLLOWED_BY_INSERTION_BASE; }
                }
            case S:
                for ( int jjj = 0; jjj < elementLength; jjj++ ) {
                    readPos++;
                }
                break;
            case D:
            case N:
                for ( int jjj = 0; jjj < elementLength; jjj++ ) {
                    alignment[alignPos] = PileupElement.DELETION_BASE;
                    alignPos++;
                }
                break;
            case M:
                for ( int jjj = 0; jjj < elementLength; jjj++ ) {
                    alignment[alignPos] = read[readPos];
                    alignPos++;
                    readPos++;
                }
                break;
            case H:
            case P:
                break;
            default:
                throw new ReviewedStingException( "Unsupported cigar operator: " + ce.getOperator() );
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
        // sam files (including those generated by maq or bwa) miss headers altogether. When
        // reading such a SAM file, reference name is set, but since there is no sequence dictionary,
        // null is always returned for referenceIndex. Let's be paranoid here, and make sure that
        // we do not call the read "unmapped" when it has only reference name set with ref. index missing
        // or vice versa.
        if ( ( r.getReferenceIndex() != null && r.getReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX
                || r.getReferenceName() != null && !r.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) )
          &&  r.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START ) return false  ;
        return true;
    }

    /**
     * Due to (unfortunate) multiple ways to indicate that read/mate is unmapped allowed by SAM format
     * specification, one may need this convenience shortcut. Checks both 'mate unmapped' flag and
     * alignment reference index/start of the mate.
     * @param r sam record for the read
     * @return true if read's mate is unmapped
     */
    public static boolean isMateUnmapped(final SAMRecord r) {
        if ( r.getMateUnmappedFlag() ) return true;

        // our life would be so much easier if all sam files followed the specs. In reality,
        // sam files (including those generated by maq or bwa) miss headers altogether. When
        // reading such a SAM file, reference name is set, but since there is no sequence dictionary,
        // null is always returned for referenceIndex. Let's be paranoid here, and make sure that
        // we do not call the read "unmapped" when it has only reference name set with ref. index missing
        // or vice versa.
        if ( ( r.getMateReferenceIndex() != null && r.getMateReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX
                || r.getMateReferenceName() != null && !r.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) )
          &&  r.getMateAlignmentStart() != SAMRecord.NO_ALIGNMENT_START ) return false  ;
        return true;
    }

    /** Returns true is read is mapped and mapped uniquely (Q>0).
     * 
     * @param read
     * @return
     */
    public static boolean isReadUniquelyMapped(SAMRecord read) {
        return ( ! AlignmentUtils.isReadUnmapped(read) ) && read.getMappingQuality() > 0;
    }

    /** Returns the array of base qualitites in the order the bases were read on the machine (i.e. always starting from
     * cycle 1). In other words, if the read is unmapped or aligned in the forward direction, the read's own base
     * qualities are returned as stored in the SAM record; if the read is aligned in the reverse direction, the array
     * of read's base qualitites is inverted (in this case new array is allocated and returned).
      * @param read
     * @return
     */
    public static byte [] getQualsInCycleOrder(SAMRecord read) {
        if ( isReadUnmapped(read) || ! read.getReadNegativeStrandFlag() ) return read.getBaseQualities();

        return Utils.reverse(read.getBaseQualities());
    }

    /** Returns the array of original base qualitites (before recalibration) in the order the bases were read on the machine (i.e. always starting from
     * cycle 1). In other words, if the read is unmapped or aligned in the forward direction, the read's own base
     * qualities are returned as stored in the SAM record; if the read is aligned in the reverse direction, the array
     * of read's base qualitites is inverted (in this case new array is allocated and returned). If no original base qualities
     * are available this method will throw a runtime exception.
      * @param read
     * @return
     */
    public static byte [] getOriginalQualsInCycleOrder(SAMRecord read) {
        if ( isReadUnmapped(read) || ! read.getReadNegativeStrandFlag() ) return read.getOriginalBaseQualities();

        return Utils.reverse(read.getOriginalBaseQualities());
    }

    /** Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
     * starting at 0-based position <code>refIndex</code> on the <code>refSeq</code> and specified by its <code>cigar</code>.
     * The last argument <code>readIndex</code> specifies 0-based position on the read where the alignment described by the
     * <code>cigar</code> starts. Usually cigars specify alignments of the whole read to the ref, so that readIndex is normally 0.
     * Use non-zero readIndex only when the alignment cigar represents alignment of a part of the read. The refIndex in this case
     * should be the position where the alignment of that part of the read starts at. In other words, both refIndex and readIndex are
     * always the positions where the cigar starts on the ref and on the read, respectively.
     *
     * If the alignment has an indel, then this method attempts moving this indel left across a stretch of repetitive bases. For instance, if the original cigar
     * specifies that (any) one AT  is deleted from a repeat sequence TATATATA, the output cigar will always mark the leftmost AT
     * as deleted. If there is no indel in the original cigar, or the indel position is determined unambiguously (i.e. inserted/deleted sequence
     * is not repeated), the original cigar is returned.
     * @param cigar structure of the original alignment
     * @param refSeq reference sequence the read is aligned to
     * @param readSeq read sequence
     * @param refIndex 0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return a cigar, in which indel is guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    public static Cigar leftAlignIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex) {

        int indexOfIndel = -1;
        for ( int i = 0; i < cigar.numCigarElements(); i++ ) {
            CigarElement ce = cigar.getCigarElement(i);
            if ( ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I ) {
                // if there is more than 1 indel, don't left align
                if ( indexOfIndel != -1 )
                    return cigar;
                indexOfIndel = i;
            }
        }

        // if there is no indel or if the alignment starts with an insertion (so that there
        // is no place on the read to move that insertion further left), we are done
        if ( indexOfIndel < 1 ) return cigar;

        final int indelLength = cigar.getCigarElement(indexOfIndel).getLength();

        byte[] altString = createIndelString(cigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);
        if ( altString == null )
            return cigar;

        Cigar newCigar = cigar;
        for ( int i = 0; i < indelLength; i++ ) {
            newCigar = moveCigarLeft(newCigar, indexOfIndel);
            byte[] newAltString = createIndelString(newCigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);

            // check to make sure we haven't run off the end of the read
            boolean reachedEndOfRead = cigarHasZeroSizeElement(newCigar);

            if ( Arrays.equals(altString, newAltString) ) {
                cigar = newCigar;
                i = -1;
                if ( reachedEndOfRead )
                    cigar = cleanUpCigar(cigar);
            }

            if ( reachedEndOfRead )
                break;
        }

        return cigar;
    }

    private static boolean cigarHasZeroSizeElement(Cigar c) {
        for ( CigarElement ce : c.getCigarElements() ) {
            if ( ce.getLength() == 0 )
                return true;
        }
        return false;
    }

    private static Cigar cleanUpCigar(Cigar c) {
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(c.numCigarElements()-1);
        for ( CigarElement ce : c.getCigarElements() ) {
            if ( ce.getLength() != 0 &&
                    (elements.size() != 0 || ce.getOperator() != CigarOperator.D) ) {               
                elements.add(ce);
            }
        }
        return new Cigar(elements);
    }

    private static Cigar moveCigarLeft(Cigar cigar, int indexOfIndel) {
        // get the first few elements
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(cigar.numCigarElements());
        for ( int i = 0; i < indexOfIndel - 1; i++)
            elements.add(cigar.getCigarElement(i));

        // get the indel element and move it left one base
        CigarElement ce = cigar.getCigarElement(indexOfIndel-1);
        elements.add(new CigarElement(ce.getLength()-1, ce.getOperator()));
        elements.add(cigar.getCigarElement(indexOfIndel));        
        if ( indexOfIndel+1 < cigar.numCigarElements() ) {
            ce = cigar.getCigarElement(indexOfIndel+1);
            elements.add(new CigarElement(ce.getLength()+1, ce.getOperator()));
        } else {
            elements.add(new CigarElement(1, CigarOperator.M));
        }

        // get the last few elements
        for ( int i = indexOfIndel + 2; i < cigar.numCigarElements(); i++)
            elements.add(cigar.getCigarElement(i));
        return new Cigar(elements);
    }

    private static byte[] createIndelString(final Cigar cigar, final int indexOfIndel, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        CigarElement indel = cigar.getCigarElement(indexOfIndel);
        int indelLength = indel.getLength();

        int totalRefBases = 0;
        for ( int i = 0; i < indexOfIndel; i++ ) {
            CigarElement ce = cigar.getCigarElement(i);
            int length = ce.getLength();
            
            switch( ce.getOperator() ) {
                case M:
                    readIndex += length;
                    refIndex += length;
                    totalRefBases += length;
                    break;
                case S:
                    readIndex += length;
                    break;
                case N:
                    refIndex += length;
                    totalRefBases += length;
                    break;
                default:
                    break;
            }
        }

        // sometimes, when there are very large known indels, we won't have enough reference sequence to cover them
        if ( totalRefBases + indelLength > refSeq.length )
            indelLength -= (totalRefBases + indelLength - refSeq.length);

        // the indel-based reference string
        byte[] alt = new byte[refSeq.length + (indelLength * (indel.getOperator() == CigarOperator.D ? -1 : 1))];

        // add the bases before the indel, making sure it's not aligned off the end of the reference
        if ( refIndex > alt.length || refIndex > refSeq.length )
            return null;
        System.arraycopy(refSeq, 0, alt, 0, refIndex);
        int currentPos = refIndex;

        // take care of the indel
        if ( indel.getOperator() == CigarOperator.D ) {
            refIndex += indelLength;
        } else {
            System.arraycopy(readSeq, readIndex, alt, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel, making sure it's not aligned off the end of the reference
        if ( refSeq.length - refIndex > alt.length - currentPos )
            return null;
        System.arraycopy(refSeq, refIndex, alt, currentPos, refSeq.length - refIndex);

        return alt;
    }
}
