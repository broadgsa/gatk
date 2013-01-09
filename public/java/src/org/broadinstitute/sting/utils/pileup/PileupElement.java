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

package org.broadinstitute.sting.utils.pileup;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 */
public class PileupElement implements Comparable<PileupElement> {
    public static final byte DELETION_BASE = BaseUtils.D;
    public static final byte DELETION_QUAL = (byte) 16;
    public static final byte A_FOLLOWED_BY_INSERTION_BASE = (byte) 87;
    public static final byte C_FOLLOWED_BY_INSERTION_BASE = (byte) 88;
    public static final byte T_FOLLOWED_BY_INSERTION_BASE = (byte) 89;
    public static final byte G_FOLLOWED_BY_INSERTION_BASE = (byte) 90;

    protected final GATKSAMRecord read;         // the read this base belongs to
    protected final int offset;                 // the offset in the bases array for this base

    private final CigarElement currentCigarElement;
    private final int currentCigarOffset;
    private final int offsetInCurrentCigar;

    /**
     * Creates a new pileup element.
     *
     * @param read              the read we are adding to the pileup
     * @param offset            the position in the read for this base. All deletions must be left aligned! (-1 is only allowed for reads starting with insertions)
     * @param isDeletion        whether or not this base is a deletion
     * @param isBeforeDeletion  whether or not this base is before a deletion
     * @param isAfterDeletion   whether or not this base is after a deletion
     * @param isBeforeInsertion whether or not this base is before an insertion
     * @param isAfterInsertion  whether or not this base is after an insertion
     * @param isNextToSoftClip  whether or not this base is next to a soft clipped base
     * @param nextEventBases    bases in event in case element comes before insertion or deletion 
     * @param nextEventLength   length of next event in case it's insertion or deletion                             
     */
    @Requires({
            "read != null",
            "offset >= -1",
            "offset <= read.getReadLength()"})
    @Deprecated
    public PileupElement(final GATKSAMRecord read, final int offset, final boolean isDeletion, final boolean isBeforeDeletion, final boolean isAfterDeletion, final boolean isBeforeInsertion, final boolean isAfterInsertion, final boolean isNextToSoftClip, final String nextEventBases, final int nextEventLength) {
        if (offset < 0 && isDeletion)
            throw new ReviewedStingException("Pileup Element cannot create a deletion with a negative offset");

        this.read = read;
        this.offset = offset;
        currentCigarElement = null;
        currentCigarOffset = offsetInCurrentCigar = -1;
    }

    public PileupElement(final GATKSAMRecord read, final int baseOffset,
                         final CigarElement currentElement, final int currentCigarOffset, final int offsetInCurrentCigar) {
        this.read = read;
        this.offset = baseOffset;
        this.currentCigarElement = currentElement;
        this.currentCigarOffset = currentCigarOffset;
        this.offsetInCurrentCigar = offsetInCurrentCigar;
    }

    /**
     * Create a new PileupElement that's a copy of toCopy
     * @param toCopy the element we want to copy
     */
    public PileupElement(final PileupElement toCopy) {
        this(toCopy.read, toCopy.offset, toCopy.currentCigarElement, toCopy.currentCigarOffset, toCopy.offsetInCurrentCigar);
    }

    @Deprecated
    public PileupElement(final GATKSAMRecord read, final int baseOffset) {
        throw new UnsupportedOperationException("please use LocusIteratorByState.createPileupForReadAndOffset instead");
    }

    public boolean isDeletion() {
        return currentCigarElement.getOperator() == CigarOperator.D;
    }

    public boolean isBeforeDeletionStart() {
        return isBeforeDeletion() && ! isDeletion();
    }

    public boolean isAfterDeletionEnd() {
        return isAfterDeletion() && ! isDeletion();
    }

    public boolean isInsertionAtBeginningOfRead() {
        return offset == -1;
    }

    @Ensures("result != null")
    public GATKSAMRecord getRead() {
        return read;
    }

    @Ensures("result == offset")
    public int getOffset() {
        return offset;
    }

    public byte getBase() {
        return getBase(offset);
    }

    public int getBaseIndex() {
        return getBaseIndex(offset);
    }

    public byte getQual() {
        return getQual(offset);
    }

    public byte getBaseInsertionQual() {
        return getBaseInsertionQual(offset);
    }

    public byte getBaseDeletionQual() {
        return getBaseDeletionQual(offset);
    }

    /**
     * Get the length of an immediately following insertion or deletion event, or 0 if no such event exists
     *
     * Only returns a positive value when this pileup element is immediately before an indel.  Being
     * immediately before a deletion means that this pileup element isn't an deletion, and that the
     * next genomic alignment for this read is a deletion.  For the insertion case, this means
     * that an insertion cigar occurs immediately after this element, between this one and the
     * next genomic position.
     *
     * Note this function may be expensive, so multiple uses should be cached by the caller
     *
     * @return length of the event (number of inserted or deleted bases), or 0
     */
    @Ensures("result >= 0")
    public int getLengthOfImmediatelyFollowingIndel() {
        final CigarElement element = getNextIndelCigarElement();
        return element == null ? 0 : element.getLength();
    }

    /**
     * Helpful function to get the immediately following cigar element, for an insertion or deletion
     *
     * if this state precedes a deletion (i.e., next position on genome) or insertion (immediately between
     * this and the next position) returns the CigarElement corresponding to this event.  Otherwise returns
     * null.
     *
     * @return a CigarElement, or null if the next alignment state ins't an insertion or deletion.
     */
    private CigarElement getNextIndelCigarElement() {
        if ( isBeforeDeletionStart() ) {
            final CigarElement element = getNextOnGenomeCigarElement();
            if ( element == null || element.getOperator() != CigarOperator.D )
                throw new IllegalStateException("Immediately before deletion but the next cigar element isn't a deletion " + element);
            return element;
        } else if ( isBeforeInsertion() ) {
            final CigarElement element = getBetweenNextPosition().get(0);
            if ( element.getOperator() != CigarOperator.I )
                throw new IllegalStateException("Immediately before insertion but the next cigar element isn't an insertion " + element);
            return element;
        } else {
            return null;
        }
    }

    /**
     * Get the bases for an insertion that immediately follows this alignment state, or null if none exists
     *
     * @see #getLengthOfImmediatelyFollowingIndel() for details on the meaning of immediately.
     *
     * If the immediately following state isn't an insertion, returns null
     *
     * @return actual sequence of inserted bases, or a null if the event is a deletion or if there is no event in the associated read.
     */
    @Ensures("result == null || result.length() == getLengthOfImmediatelyFollowingIndel()")
    public String getBasesOfImmediatelyFollowingInsertion() {
        final CigarElement element = getNextIndelCigarElement();
        if ( element != null && element.getOperator() == CigarOperator.I ) {
            final int getFrom = offset + 1;
            final byte[] bases = Arrays.copyOfRange(read.getReadBases(), getFrom, getFrom + element.getLength());
            return new String(bases);
        } else
            return null;
    }

    public int getMappingQual() {
        return read.getMappingQuality();
    }

    @Ensures("result != null")
    public String toString() {
        return String.format("%s @ %d = %c Q%d", getRead().getReadName(), getOffset(), (char) getBase(), getQual());
    }

    protected byte getBase(final int offset) {
        return (isDeletion() || isInsertionAtBeginningOfRead()) ? DELETION_BASE : read.getReadBases()[offset];
    }

    protected int getBaseIndex(final int offset) {
        return BaseUtils.simpleBaseToBaseIndex((isDeletion() || isInsertionAtBeginningOfRead()) ? DELETION_BASE : read.getReadBases()[offset]);
    }

    protected byte getQual(final int offset) {
        return (isDeletion() || isInsertionAtBeginningOfRead()) ? DELETION_QUAL : read.getBaseQualities()[offset];
    }

    protected byte getBaseInsertionQual(final int offset) {
        return (isDeletion() || isInsertionAtBeginningOfRead()) ? DELETION_QUAL : read.getBaseInsertionQualities()[offset];
    }

    protected byte getBaseDeletionQual(final int offset) {
        return (isDeletion() || isInsertionAtBeginningOfRead()) ? DELETION_QUAL : read.getBaseDeletionQualities()[offset];
    }

    @Override
    public int compareTo(final PileupElement pileupElement) {
        if (offset < pileupElement.offset)
            return -1;
        else if (offset > pileupElement.offset)
            return 1;
        else if (read.getAlignmentStart() < pileupElement.read.getAlignmentStart())
            return -1;
        else if (read.getAlignmentStart() > pileupElement.read.getAlignmentStart())
            return 1;
        else
            return 0;
    }

    // --------------------------------------------------------------------------
    //
    // Reduced read accessors
    //
    // --------------------------------------------------------------------------

    /**
     * Returns the number of elements in the pileup element.
     *
     * Unless this is a reduced read, the number of elements in a pileup element is one. In the event of
     * this being a reduced read and a deletion, we return the average number of elements between the left
     * and right elements to the deletion. We assume the deletion to be left aligned.
     *
     * @return
     */
    public int getRepresentativeCount() {
        int representativeCount = 1;

        if (read.isReducedRead() && !isInsertionAtBeginningOfRead())     {
            if (isDeletion() && (offset + 1 >= read.getReadLength()) )  // deletion in the end of the read
                throw new UserException.MalformedBAM(read, String.format("Adjacent I/D events in read %s -- cigar: %s", read.getReadName(), read.getCigarString()));

            representativeCount = (isDeletion()) ? MathUtils.fastRound((read.getReducedCount(offset) + read.getReducedCount(offset + 1)) / 2.0) : read.getReducedCount(offset);
        }
        return representativeCount;
    }

    public CigarElement getCurrentCigarElement() {
        return currentCigarElement;
    }

    public int getCurrentCigarOffset() {
        return currentCigarOffset;
    }

    public int getOffsetInCurrentCigar() {
        return offsetInCurrentCigar;
    }

    public LinkedList<CigarElement> getBetweenPrevPosition() {
        return atStartOfCurrentCigar() ? getBetween(-1) : EMPTY_LINKED_LIST;
    }

    public LinkedList<CigarElement> getBetweenNextPosition() {
        return atEndOfCurrentCigar() ? getBetween(1) : EMPTY_LINKED_LIST;
    }

    // TODO -- can I make this unmodifable?
    private final static LinkedList<CigarElement> EMPTY_LINKED_LIST = new LinkedList<CigarElement>();

    private final static EnumSet<CigarOperator> ON_GENOME_OPERATORS =
            EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X, CigarOperator.D);

    private LinkedList<CigarElement> getBetween(final int increment) {
        LinkedList<CigarElement> elements = null;
        final int nCigarElements = read.getCigarLength();
        for ( int i = currentCigarOffset + increment; i >= 0 && i < nCigarElements; i += increment) {
            final CigarElement elt = read.getCigar().getCigarElement(i);
            if ( ON_GENOME_OPERATORS.contains(elt.getOperator()) )
                break;
            else {
                // optimization: don't allocate list if not necessary
                if ( elements == null )
                    elements = new LinkedList<CigarElement>();

                if ( increment > 0 )
                    // to keep the list in the right order, if we are incrementing positively add to the end
                    elements.add(elt);
                else
                    // counting down => add to front
                    elements.addFirst(elt);
            }
        }

        // optimization: elements is null because nothing got added, just return the empty list
        return elements == null ? EMPTY_LINKED_LIST : elements;
    }

    public CigarElement getPreviousOnGenomeCigarElement() {
        return getNeighboringOnGenomeCigarElement(-1);
    }

    public CigarElement getNextOnGenomeCigarElement() {
        return getNeighboringOnGenomeCigarElement(1);
    }

    private CigarElement getNeighboringOnGenomeCigarElement(final int increment) {
        final int nCigarElements = read.getCigarLength();

        for ( int i = currentCigarOffset + increment; i >= 0 && i < nCigarElements; i += increment) {
            final CigarElement elt = read.getCigar().getCigarElement(i);
            if ( ON_GENOME_OPERATORS.contains(elt.getOperator()) )
                return elt;
        }

        // getting here means that you didn't find anything
        return null;
    }

    private boolean hasOperator(final CigarElement maybeCigarElement, final CigarOperator toMatch) {
        return maybeCigarElement != null && maybeCigarElement.getOperator() == toMatch;
    }

    public boolean isAfterDeletion() { return atStartOfCurrentCigar() && hasOperator(getPreviousOnGenomeCigarElement(), CigarOperator.D); }
    public boolean isBeforeDeletion() { return atEndOfCurrentCigar() && hasOperator(getNextOnGenomeCigarElement(), CigarOperator.D); }
    public boolean isAfterInsertion() { return isAfter(getBetweenPrevPosition(), CigarOperator.I); }
    public boolean isBeforeInsertion() { return isBefore(getBetweenNextPosition(), CigarOperator.I); }

    public boolean isAfterSoftClip() { return isAfter(getBetweenPrevPosition(), CigarOperator.S); }
    public boolean isBeforeSoftClip() { return isBefore(getBetweenNextPosition(), CigarOperator.S); }
    public boolean isNextToSoftClip() { return isAfterSoftClip() || isBeforeSoftClip(); }

    public boolean atEndOfCurrentCigar() {
        return offsetInCurrentCigar == currentCigarElement.getLength() - 1;
    }

    public boolean atStartOfCurrentCigar() {
        return offsetInCurrentCigar == 0;
    }

    private boolean isAfter(final LinkedList<CigarElement> elements, final CigarOperator op) {
        return ! elements.isEmpty() && elements.peekLast().getOperator() == op;
    }

    private boolean isBefore(final List<CigarElement> elements, final CigarOperator op) {
        return ! elements.isEmpty() && elements.get(0).getOperator() == op;
    }
}