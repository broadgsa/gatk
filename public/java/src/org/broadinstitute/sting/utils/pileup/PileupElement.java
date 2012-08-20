package org.broadinstitute.sting.utils.pileup;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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
    protected final boolean isDeletion;         // is this base a deletion
    protected final boolean isBeforeDeletedBase;   // is the base to the right of this base an deletion
    protected final boolean isAfterDeletedBase;    // is the base to the left of this base a deletion
    protected final boolean isBeforeInsertion;  // is the base to the right of this base an insertion
    protected final boolean isAfterInsertion;   // is the base to the left of this base an insertion
    protected final boolean isNextToSoftClip;   // is this base either before or after a soft clipped base
    protected final int eventLength;            // what is the length of the event (insertion or deletion) *after* this base
    protected final String eventBases;          // if it is a deletion, we do not have information about the actual deleted bases in the read itself, so we fill the string with D's; for insertions we keep actual inserted bases

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
    public PileupElement(final GATKSAMRecord read, final int offset, final boolean isDeletion, final boolean isBeforeDeletion, final boolean isAfterDeletion, final boolean isBeforeInsertion, final boolean isAfterInsertion, final boolean isNextToSoftClip, final String nextEventBases, final int nextEventLength) {
        if (offset < 0 && isDeletion)
            throw new ReviewedStingException("Pileup Element cannot create a deletion with a negative offset");

        this.read = read;
        this.offset = offset;
        this.isDeletion = isDeletion;
        this.isBeforeDeletedBase = isBeforeDeletion;
        this.isAfterDeletedBase = isAfterDeletion;
        this.isBeforeInsertion = isBeforeInsertion;
        this.isAfterInsertion = isAfterInsertion;
        this.isNextToSoftClip = isNextToSoftClip;
        if (isBeforeInsertion)
            eventBases = nextEventBases;
        else
            eventBases = null;                  // ignore argument in any other case
        if (isBeforeDeletion || isBeforeInsertion)
            eventLength = nextEventLength;
        else
            eventLength = -1;
    }

    public PileupElement(final GATKSAMRecord read, final int offset, final boolean isDeletion, final boolean isBeforeDeletion, final boolean isAfterDeletion, final boolean isBeforeInsertion, final boolean isAfterInsertion, final boolean isNextToSoftClip) {
        this(read, offset, isDeletion, isBeforeDeletion, isAfterDeletion, isBeforeInsertion, isAfterInsertion, isNextToSoftClip, null, -1);
    }
    public boolean isDeletion() {
        return isDeletion;
    }

    public boolean isBeforeDeletedBase() {
        return isBeforeDeletedBase;
    }

    public boolean isAfterDeletedBase() {
        return isAfterDeletedBase;
    }

    public boolean isBeforeDeletionStart() {
        return isBeforeDeletedBase && !isDeletion;
    }

    public boolean isAfterDeletionEnd() {
        return isAfterDeletedBase && !isDeletion;
    }

    public boolean isBeforeInsertion() {
        return isBeforeInsertion;
    }

    public boolean isAfterInsertion() {
        return isAfterInsertion;
    }

    public boolean isNextToSoftClip() {
        return isNextToSoftClip;
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
     * @return length of the event (number of inserted or deleted bases
     */
    public int getEventLength() {
        return eventLength;
    }

    /**
     * @return actual sequence of inserted bases, or a null if the event is a deletion or if there is no event in the associated read.
     */
    public String getEventBases() {
        return eventBases;
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

}