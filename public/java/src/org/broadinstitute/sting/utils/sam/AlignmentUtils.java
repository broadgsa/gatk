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

package org.broadinstitute.sting.utils.sam;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.recalibration.EventType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;


public final class AlignmentUtils {
    private final static EnumSet<CigarOperator> ALIGNED_TO_GENOME_OPERATORS = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X);
    private final static EnumSet<CigarOperator> ALIGNED_TO_GENOME_PLUS_SOFTCLIPS = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X, CigarOperator.S);

    // cannot be instantiated
    private AlignmentUtils() { }

    public static class MismatchCount {
        public int numMismatches = 0;
        public long mismatchQualities = 0;
    }

    public static long mismatchingQualities(GATKSAMRecord r, byte[] refSeq, int refIndex) {
        return getMismatchCount(r, refSeq, refIndex).mismatchQualities;
    }

    public static MismatchCount getMismatchCount(GATKSAMRecord r, byte[] refSeq, int refIndex) {
        return getMismatchCount(r, refSeq, refIndex, 0, r.getReadLength());
    }

    // todo -- this code and mismatchesInRefWindow should be combined and optimized into a single
    // todo -- high performance implementation.  We can do a lot better than this right now

    /**
     * Count how many bases mismatch the reference.  Indels are not considered mismatching.
     *
     * @param r                   the sam record to check against
     * @param refSeq              the byte array representing the reference sequence
     * @param refIndex            the index in the reference byte array of the read's first base (the reference index is matching the alignment start, there may be tons of soft-clipped bases before/after that so it's wrong to compare with getReadLength() here.)
     * @param startOnRead         the index in the read's bases from which we start counting
     * @param nReadBases          the number of bases after (but including) startOnRead that we check
     * @return non-null object representing the mismatch count
     */
    @Ensures("result != null")
    public static MismatchCount getMismatchCount(GATKSAMRecord r, byte[] refSeq, int refIndex, int startOnRead, int nReadBases) {
        if ( r == null ) throw new IllegalArgumentException("attempting to calculate the mismatch count from a read that is null");
        if ( refSeq == null ) throw new IllegalArgumentException("attempting to calculate the mismatch count with a reference sequence that is null");
        if ( refIndex < 0 ) throw new IllegalArgumentException("attempting to calculate the mismatch count with a reference index that is negative");
        if ( startOnRead < 0 ) throw new IllegalArgumentException("attempting to calculate the mismatch count with a read start that is negative");
        if ( nReadBases < 0 ) throw new IllegalArgumentException("attempting to calculate the mismatch count for a negative number of read bases");
        if ( refSeq.length - refIndex < (r.getAlignmentEnd() - r.getAlignmentStart()) )
            throw new IllegalArgumentException("attempting to calculate the mismatch count against a reference string that is smaller than the read");

        MismatchCount mc = new MismatchCount();

        int readIdx = 0;
        final int endOnRead = startOnRead + nReadBases - 1; // index of the last base on read we want to count (note we are including soft-clipped bases with this math)
        final byte[] readSeq = r.getReadBases();
        final Cigar c = r.getCigar();
        final byte[] readQuals = r.getBaseQualities();
        for (final CigarElement ce : c.getCigarElements()) {

            if (readIdx > endOnRead)
                break;

            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case X:
                    mc.numMismatches += elementLength;
                    for (int j = 0; j < elementLength; j++)
                        mc.mismatchQualities += readQuals[readIdx+j];
                case EQ:
                    refIndex += elementLength;
                    readIdx += elementLength;
                break;
                case M:
                    for (int j = 0; j < elementLength; j++, refIndex++, readIdx++) {
                        if (refIndex >= refSeq.length)
                            continue;                      // TODO : It should never happen, we should throw exception here
                        if (readIdx < startOnRead) continue;
                        if (readIdx > endOnRead) break;
                        byte refChr = refSeq[refIndex];
                        byte readChr = readSeq[readIdx];
                        // Note: we need to count X/N's as mismatches because that's what SAM requires
                        //if ( BaseUtils.simpleBaseToBaseIndex(readChr) == -1 ||
                        //     BaseUtils.simpleBaseToBaseIndex(refChr)  == -1 )
                        //    continue; // do not count Ns/Xs/etc ?
                        if (readChr != refChr) {
                            mc.numMismatches++;
                            mc.mismatchQualities += readQuals[readIdx];
                        }
                    }
                    break;
                case I:
                case S:
                    readIdx += elementLength;
                    break;
                case D:
                case N:
                    refIndex += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new ReviewedStingException("The " + ce.getOperator() + " cigar element is not currently supported");
            }

        }
        return mc;
    }

    /**
     * Returns number of alignment blocks (continuous stretches of aligned bases) in the specified alignment.
     * This method follows closely the SAMRecord::getAlignmentBlocks() implemented in samtools library, but
     * it only counts blocks without actually allocating and filling the list of blocks themselves. Hence, this method is
     * a much more efficient alternative to r.getAlignmentBlocks.size() in the situations when this number is all that is needed.
     * Formally, this method simply returns the number of M elements in the cigar.
     *
     * @param r alignment
     * @return number of continuous alignment blocks (i.e. 'M' elements of the cigar; all indel and clipping elements are ignored).
     */
    @Ensures("result >= 0")
    public static int getNumAlignmentBlocks(final SAMRecord r) {
        if ( r == null ) throw new IllegalArgumentException("read cannot be null");
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;

        int n = 0;
        for (final CigarElement e : cigar.getCigarElements()) {
            if (ALIGNED_TO_GENOME_OPERATORS.contains(e.getOperator()))
                n++;
        }

        return n;
    }


    /**
     * Get the number of bases aligned to the genome, including soft clips
     *
     * If read is not mapped (i.e., doesn't have a cigar) returns 0
     *
     * @param r a non-null GATKSAMRecord
     * @return the number of bases aligned to the genome in R, including soft clipped bases
     */
    public static int getNumAlignedBasesCountingSoftClips(final GATKSAMRecord r) {
        int n = 0;
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;

        for (final CigarElement e : cigar.getCigarElements())
            if (ALIGNED_TO_GENOME_PLUS_SOFTCLIPS.contains(e.getOperator()))
                n += e.getLength();

        return n;
    }

    /**
     * Count the number of bases hard clipped from read
     *
     * If read's cigar is null, return 0
     *
     * @param r a non-null read
     * @return a positive integer
     */
    @Ensures("result >= 0")
    public static int getNumHardClippedBases(final SAMRecord r) {
        if ( r == null ) throw new IllegalArgumentException("Read cannot be null");

        int n = 0;
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;

        for (final CigarElement e : cigar.getCigarElements())
            if (e.getOperator() == CigarOperator.H)
                n += e.getLength();

        return n;
    }

    /**
     * Calculate the number of bases that are soft clipped in read with quality score greater than threshold
     *
     * Handles the case where the cigar is null (i.e., the read is unmapped), returning 0
     *
     * @param read a non-null GATKSAMRecord.
     * @param qualThreshold consider bases with quals > this value as high quality.  Must be >= 0
     * @return positive integer
     */
    @Ensures("result >= 0")
    public static int calcNumHighQualitySoftClips( final GATKSAMRecord read, final byte qualThreshold ) {
        if ( read == null ) throw new IllegalArgumentException("Read cannot be null");
        if ( qualThreshold < 0 ) throw new IllegalArgumentException("Expected qualThreshold to be a positive byte but saw " + qualThreshold);

        if ( read.getCigar() == null ) // the read is unmapped
            return 0;

        final byte[] qual = read.getBaseQualities( EventType.BASE_SUBSTITUTION );

        int numHQSoftClips = 0;
        int alignPos = 0;
        for ( final CigarElement ce : read.getCigar().getCigarElements() ) {
            final int elementLength = ce.getLength();

            switch( ce.getOperator() ) {
                case S:
                    for( int jjj = 0; jjj < elementLength; jjj++ ) {
                        if( qual[alignPos++] > qualThreshold ) { numHQSoftClips++; }
                    }
                    break;
                case M: case I: case EQ: case X:
                    alignPos += elementLength;
                    break;
                case H: case P: case D: case N:
                    break;
                default:
                    throw new IllegalStateException("Unsupported cigar operator: " + ce.getOperator());
            }
        }

        return numHQSoftClips;
    }

    public static int calcAlignmentByteArrayOffset(final Cigar cigar, final PileupElement pileupElement, final int alignmentStart, final int refLocus) {
        return calcAlignmentByteArrayOffset( cigar, pileupElement.getOffset(), pileupElement.isDeletion(), alignmentStart, refLocus );
    }

    /**
     * Calculate the index into the read's bases of the beginning of the encompassing cigar element for a given cigar and offset
     *
     * @param cigar            the read's CIGAR -- cannot be null
     * @param offset           the offset to use for the calculation or -1 if in the middle of a deletion
     * @param isDeletion       are we in the middle of a deletion?
     * @param alignmentStart   the alignment start of the read
     * @param refLocus         the reference position of the offset
     * @return a non-negative int index
     */
    @Ensures("result >= 0")
    public static int calcAlignmentByteArrayOffset(final Cigar cigar, final int offset, final boolean isDeletion, final int alignmentStart, final int refLocus) {
        if ( cigar == null ) throw new IllegalArgumentException("attempting to find the alignment position from a CIGAR that is null");
        if ( offset < -1 ) throw new IllegalArgumentException("attempting to find the alignment position with an offset that is negative (and not -1)");
        if ( alignmentStart < 0 ) throw new IllegalArgumentException("attempting to find the alignment position from an alignment start that is negative");
        if ( refLocus < 0 ) throw new IllegalArgumentException("attempting to find the alignment position from a reference position that is negative");
        if ( offset >= cigar.getReadLength() ) throw new IllegalArgumentException("attempting to find the alignment position of an offset than is larger than the read length");

        int pileupOffset = offset;

        // Reassign the offset if we are in the middle of a deletion because of the modified representation of the read bases
        if (isDeletion) {
            pileupOffset = refLocus - alignmentStart;
            final CigarElement ce = cigar.getCigarElement(0);
            if (ce.getOperator() == CigarOperator.S) {
                pileupOffset += ce.getLength();
            }
        }

        int pos = 0;
        int alignmentPos = 0;

        for (int iii = 0; iii < cigar.numCigarElements(); iii++) {
            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch (ce.getOperator()) {
                case I:
                case S: // TODO -- I don't think that soft clips should be treated the same as inserted bases here. Investigation needed.
                    pos += elementLength;
                    if (pos >= pileupOffset) {
                        return alignmentPos;
                    }
                    break;
                case D:
                    if (!isDeletion) {
                        alignmentPos += elementLength;
                    } else {
                        if (pos + elementLength - 1 >= pileupOffset) {
                            return alignmentPos + (pileupOffset - pos);
                        } else {
                            pos += elementLength;
                            alignmentPos += elementLength;
                        }
                    }
                    break;
                case M:
                case EQ:
                case X:
                    if (pos + elementLength - 1 >= pileupOffset) {
                        return alignmentPos + (pileupOffset - pos);
                    } else {
                        pos += elementLength;
                        alignmentPos += elementLength;
                    }
                    break;
                case H:
                case P:
                case N:
                    break;
                default:
                    throw new ReviewedStingException("Unsupported cigar operator: " + ce.getOperator());
            }
        }

        return alignmentPos;
    }

    /**
     * Generate an array of bases for just those that are aligned to the reference (i.e. no clips or insertions)
     *
     * @param cigar            the read's CIGAR -- cannot be null
     * @param read             the read's base array
     * @return a non-null array of bases (bytes)
     */
    @Ensures("result != null")
    public static byte[] readToAlignmentByteArray(final Cigar cigar, final byte[] read) {
        if ( cigar == null ) throw new IllegalArgumentException("attempting to generate an alignment from a CIGAR that is null");
        if ( read == null ) throw new IllegalArgumentException("attempting to generate an alignment from a read sequence that is null");

        final int alignmentLength = cigar.getReferenceLength();
        final byte[] alignment = new byte[alignmentLength];
        int alignPos = 0;
        int readPos = 0;
        for (int iii = 0; iii < cigar.numCigarElements(); iii++) {

            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch (ce.getOperator()) {
                case I:
                    if (alignPos > 0) {
                        final int prevPos = alignPos - 1;
                        if (alignment[prevPos] == BaseUtils.Base.A.base) {
                            alignment[prevPos] = PileupElement.A_FOLLOWED_BY_INSERTION_BASE;
                        } else if (alignment[prevPos] == BaseUtils.Base.C.base) {
                            alignment[prevPos] = PileupElement.C_FOLLOWED_BY_INSERTION_BASE;
                        } else if (alignment[prevPos] == BaseUtils.Base.T.base) {
                            alignment[prevPos] = PileupElement.T_FOLLOWED_BY_INSERTION_BASE;
                        } else if (alignment[prevPos] == BaseUtils.Base.G.base) {
                            alignment[prevPos] = PileupElement.G_FOLLOWED_BY_INSERTION_BASE;
                        }
                    }
                case S:
                    readPos += elementLength;
                    break;
                case D:
                case N:
                    for (int jjj = 0; jjj < elementLength; jjj++) {
                        alignment[alignPos++] = PileupElement.DELETION_BASE;
                    }
                    break;
                case M:
                case EQ:
                case X:
                    for (int jjj = 0; jjj < elementLength; jjj++) {
                        alignment[alignPos++] = read[readPos++];
                    }
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new ReviewedStingException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return alignment;
    }

    /**
     * Returns true if the read does not belong to a contig, i.e. it's location is GenomeLoc.UNMAPPED.
     * NOTE: A read can have a mapped GenomeLoc and still have an unmapped flag!
     *
     * @param r record
     * @return true if read is unmapped to a genome loc
     */
    public static boolean isReadGenomeLocUnmapped(final SAMRecord r) {
        return SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(r.getReferenceName());
    }

    /**
     * Due to (unfortunate) multiple ways to indicate that read is unmapped allowed by SAM format
     * specification, one may need this convenience shortcut. Checks both 'read unmapped' flag and
     * alignment reference index/start.
     *
     * Our life would be so much easier if all sam files followed the specs. In reality,
     * sam files (including those generated by maq or bwa) miss headers altogether. When
     * reading such a SAM file, reference name is set, but since there is no sequence dictionary,
     * null is always returned for referenceIndex. Let's be paranoid here, and make sure that
     * we do not call the read "unmapped" when it has only reference name set with ref. index missing
     * or vice versa.
     *
     * @param r a non-null record
     * @return true if read is unmapped
     */
    public static boolean isReadUnmapped(final SAMRecord r) {
        if ( r == null )
            throw new IllegalArgumentException("Read cannot be null");

        return r.getReadUnmappedFlag() ||
               !((r.getReferenceIndex() != null && r.getReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX ||
                  r.getReferenceName() != null && !r.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) &&
                 r.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START);

    }

    /**
     * Need a well-formed, consolidated Cigar string so that the left aligning code works properly.
     * For example, 1M1M1M1D2M1M --> 3M1D3M
     * If the given cigar is empty then the returned cigar will also be empty
     * @param c the cigar to consolidate
     * @return  a non-null cigar with consecutive matching operators merged into single operators.
     */
    @Ensures({"result != null"})
    public static Cigar consolidateCigar( final Cigar c ) {
        if( c == null ) { throw new IllegalArgumentException("Cigar cannot be null"); }
        if( c.isEmpty() ) { return c; }

        final Cigar returnCigar = new Cigar();
        int sumLength = 0;
        for( int iii = 0; iii < c.numCigarElements(); iii++ ) {
            sumLength += c.getCigarElement(iii).getLength();
            if( iii == c.numCigarElements() - 1 || !c.getCigarElement(iii).getOperator().equals(c.getCigarElement(iii+1).getOperator())) { // at the end so finish the current element
                returnCigar.add(new CigarElement(sumLength, c.getCigarElement(iii).getOperator()));
                sumLength = 0;
            }
        }
        return returnCigar;
    }

    /**
     * Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
     * starting at 0-based position <code>refIndex</code> on the <code>refSeq</code> and specified by its <code>cigar</code>.
     * The last argument <code>readIndex</code> specifies 0-based position on the read where the alignment described by the
     * <code>cigar</code> starts. Usually cigars specify alignments of the whole read to the ref, so that readIndex is normally 0.
     * Use non-zero readIndex only when the alignment cigar represents alignment of a part of the read. The refIndex in this case
     * should be the position where the alignment of that part of the read starts at. In other words, both refIndex and readIndex are
     * always the positions where the cigar starts on the ref and on the read, respectively.
     * <p/>
     * If the alignment has one or more indels, this method attempts to move them left across a stretch of repetitive bases.
     * For instance, if the original cigar specifies that (any) one AT is deleted from a repeat sequence TATATATA, the output
     * cigar will always mark the leftmost AT as deleted. If there is no indel in the original cigar or if the indel position
     * is determined unambiguously (i.e. inserted/deleted sequence is not repeated), the original cigar is returned.
     *
     * Note that currently we do not actually support the case where there is more than one indel in the alignment.  We will throw
     * an exception if there is -- unless the
     *
     * @param cigar     structure of the original alignment
     * @param refSeq    reference sequence the read is aligned to
     * @param readSeq   read sequence
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @param doNotThrowExceptionForMultipleIndels  if true we will not throw an exception if we encounter multiple indels in the alignment will instead will return the original cigar
     * @return a non-null cigar, in which the indels are guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    @Ensures("result != null")
    public static Cigar leftAlignIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex, final boolean doNotThrowExceptionForMultipleIndels) {
        ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

        final int numIndels = countIndelElements(cigar);
        if ( numIndels == 0 )
            return cigar;
        if ( numIndels == 1 )
            return leftAlignSingleIndel(cigar, refSeq, readSeq, refIndex, readIndex);

        // if we got here then there is more than 1 indel in the alignment
        if ( doNotThrowExceptionForMultipleIndels )
            return cigar;

        throw new UnsupportedOperationException("attempting to left align a CIGAR that has more than 1 indel in its alignment but this functionality has not been implemented yet");
    }

    private static void ensureLeftAlignmentHasGoodArguments(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex) {
        if ( cigar == null ) throw new IllegalArgumentException("attempting to left align a CIGAR that is null");
        if ( refSeq == null ) throw new IllegalArgumentException("attempting to left align a reference sequence that is null");
        if ( readSeq == null ) throw new IllegalArgumentException("attempting to left align a read sequence that is null");
        if ( refIndex < 0 ) throw new IllegalArgumentException("attempting to left align with a reference index less than 0");
        if ( readIndex < 0 ) throw new IllegalArgumentException("attempting to left align with a read index less than 0");
    }

    /**
     * Counts the number of I/D operators
     *
     * @param cigar   cigar to check -- cannot be null
     * @return  non-negative count of indel operators
     */
    @Requires("cigar != null")
    @Ensures("result >= 0")
    private static int countIndelElements(final Cigar cigar) {
        int indelCount = 0;
        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I )
                indelCount++;
        }
        return indelCount;
    }

    /**
     * See the documentation for AlignmentUtils.leftAlignIndel() for more details.
     *
     * This flavor of the left alignment works if and only if the alignment has one - and only one - indel.
     * An exception is thrown if there are no indels or more than 1 indel in the alignment.
     *
     * @param cigar     structure of the original alignment -- cannot be null
     * @param refSeq    reference sequence the read is aligned to
     * @param readSeq   read sequence
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return a non-null cigar, in which the single indel is guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    @Ensures("result != null")
    public static Cigar leftAlignSingleIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex) {
        ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

        int indexOfIndel = -1;
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                // if there is more than 1 indel, exception out
                if (indexOfIndel != -1)
                    throw new IllegalArgumentException("attempting to left align a CIGAR that has more than 1 indel in its alignment");
                indexOfIndel = i;
            }
        }

        // if there is no indel, exception out
        if ( indexOfIndel == -1 )
            throw new IllegalArgumentException("attempting to left align a CIGAR that has no indels in its alignment");
        // if the alignment starts with an insertion (so that there is no place on the read to move that insertion further left), we are done
        if ( indexOfIndel == 0 )
            return cigar;

        final int indelLength = cigar.getCigarElement(indexOfIndel).getLength();

        byte[] altString = createIndelString(cigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);
        if (altString == null)
            return cigar;

        Cigar newCigar = cigar;
        for (int i = 0; i < indelLength; i++) {
            newCigar = moveCigarLeft(newCigar, indexOfIndel);
            byte[] newAltString = createIndelString(newCigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);

            // check to make sure we haven't run off the end of the read
            boolean reachedEndOfRead = cigarHasZeroSizeElement(newCigar);

            if (Arrays.equals(altString, newAltString)) {
                cigar = newCigar;
                i = -1;
                if (reachedEndOfRead)
                    cigar = cleanUpCigar(cigar);
            }

            if (reachedEndOfRead)
                break;
        }

        return cigar;
    }

    /**
     * Does one of the elements in cigar have a 0 length?
     *
     * @param c a non-null cigar
     * @return true if any element has 0 size
     */
    @Requires("c != null")
    protected static boolean cigarHasZeroSizeElement(final Cigar c) {
        for (final CigarElement ce : c.getCigarElements()) {
            if (ce.getLength() == 0)
                return true;
        }
        return false;
    }

    /**
     * Clean up the incoming cigar
     *
     * Removes elements with zero size
     * Clips away beginning deletion operators
     *
     * @param c the cigar string we want to clean up
     * @return a newly allocated, cleaned up Cigar
     */
    @Requires("c != null")
    @Ensures("result != null")
    private static Cigar cleanUpCigar(final Cigar c) {
        final List<CigarElement> elements = new ArrayList<CigarElement>(c.numCigarElements() - 1);

        for (final CigarElement ce : c.getCigarElements()) {
            if (ce.getLength() != 0 && (! elements.isEmpty() || ce.getOperator() != CigarOperator.D)) {
                elements.add(ce);
            }
        }

        return new Cigar(elements);
    }

    /**
     * Move the indel in a given cigar string one base to the left
     *
     * @param cigar          original cigar
     * @param indexOfIndel   the index of the indel cigar element
     * @return non-null cigar with indel moved one base to the left
     */
    @Requires("cigar != null && indexOfIndel >= 0 && indexOfIndel < cigar.numCigarElements()")
    @Ensures("result != null")
    private static Cigar moveCigarLeft(Cigar cigar, int indexOfIndel) {
        // get the first few elements
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(cigar.numCigarElements());
        for (int i = 0; i < indexOfIndel - 1; i++)
            elements.add(cigar.getCigarElement(i));

        // get the indel element and move it left one base
        CigarElement ce = cigar.getCigarElement(indexOfIndel - 1);
        elements.add(new CigarElement(Math.max(ce.getLength() - 1, 0), ce.getOperator()));
        elements.add(cigar.getCigarElement(indexOfIndel));
        if (indexOfIndel + 1 < cigar.numCigarElements()) {
            ce = cigar.getCigarElement(indexOfIndel + 1);
            elements.add(new CigarElement(ce.getLength() + 1, ce.getOperator()));
        } else {
            elements.add(new CigarElement(1, CigarOperator.M));
        }

        // get the last few elements
        for (int i = indexOfIndel + 2; i < cigar.numCigarElements(); i++)
            elements.add(cigar.getCigarElement(i));
        return new Cigar(elements);
    }

    /**
     * Create the string (really a byte array) representation of an indel-containing cigar against the reference.
     *
     * @param cigar             the indel-containing cigar
     * @param indexOfIndel      the index of the indel cigar element
     * @param refSeq            the reference sequence
     * @param readSeq           the read sequence for the cigar
     * @param refIndex          the starting reference index into refSeq
     * @param readIndex         the starting read index into readSeq
     * @return non-null byte array which is the indel representation against the reference
     */
    @Requires("cigar != null && indexOfIndel >= 0 && indexOfIndel < cigar.numCigarElements() && refSeq != null && readSeq != null && refIndex >= 0 && readIndex >= 0")
    @Ensures("result != null")
    private static byte[] createIndelString(final Cigar cigar, final int indexOfIndel, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        CigarElement indel = cigar.getCigarElement(indexOfIndel);
        int indelLength = indel.getLength();

        int totalRefBases = 0;
        for (int i = 0; i < indexOfIndel; i++) {
            CigarElement ce = cigar.getCigarElement(i);
            int length = ce.getLength();

            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
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
        if (totalRefBases + indelLength > refSeq.length)
            indelLength -= (totalRefBases + indelLength - refSeq.length);

        // the indel-based reference string
        byte[] alt = new byte[refSeq.length + (indelLength * (indel.getOperator() == CigarOperator.D ? -1 : 1))];

        // add the bases before the indel, making sure it's not aligned off the end of the reference
        if (refIndex > alt.length || refIndex > refSeq.length)
            return null;
        System.arraycopy(refSeq, 0, alt, 0, refIndex);
        int currentPos = refIndex;

        // take care of the indel
        if (indel.getOperator() == CigarOperator.D) {
            refIndex += indelLength;
        } else {
            System.arraycopy(readSeq, readIndex, alt, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel, making sure it's not aligned off the end of the reference
        if (refSeq.length - refIndex > alt.length - currentPos)
            return null;
        System.arraycopy(refSeq, refIndex, alt, currentPos, refSeq.length - refIndex);

        return alt;
    }
}
