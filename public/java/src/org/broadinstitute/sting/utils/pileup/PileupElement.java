package org.broadinstitute.sting.utils.pileup;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
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

    protected final GATKSAMRecord read;
    protected final int offset;
    protected final boolean isDeletion;
    protected final boolean isBeforeInsertion;
    protected final boolean isSoftClipped;


    /**
     * Creates a new pileup element.
     *
     * @param read              the read we are adding to the pileup
     * @param offset            the position in the read for this base. All deletions must be left aligned! (-1 is only allowed for reads starting with insertions)
     * @param isDeletion        whether or not this base is a deletion
     * @param isBeforeInsertion whether or not this base is before an insertion
     * @param isSoftClipped     whether or not this base was softclipped
     */
    @Requires({
            "read != null",
            "offset >= -1",
            "offset <= read.getReadLength()"})
    public PileupElement(final GATKSAMRecord read, final int offset, final boolean isDeletion, final boolean isBeforeInsertion, final boolean isSoftClipped) {
        if (offset < 0 && isDeletion)
            throw new ReviewedStingException("Pileup Element cannot create a deletion with a negative offset");

        this.read = read;
        this.offset = offset;
        this.isDeletion = isDeletion;
        this.isBeforeInsertion = isBeforeInsertion;
        this.isSoftClipped = isSoftClipped;
    }

    public boolean isDeletion() {
        return isDeletion;
    }

    public boolean isBeforeInsertion() {
        return isBeforeInsertion;
    }

    public boolean isSoftClipped() {
        return isSoftClipped;
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

//    public boolean isReducedRead() {
//        return read.isReducedRead();
//    }

    /**
     * Returns the number of elements in the pileup element.
     * <p/>
     * Unless this is a reduced read, the number of elements in a pileup element is one. In the event of
     * this being a reduced read and a deletion, we return the average number of elements between the left
     * and right elements to the deletion. We assume the deletion to be left aligned.
     *
     * @return
     */
    public int getRepresentativeCount() {
        int representativeCount = 1;

        if (read.isReducedRead() && !isInsertionAtBeginningOfRead())
            representativeCount = (isDeletion()) ? Math.round((read.getReducedCount(offset) + read.getReducedCount(offset + 1)) / 2) : read.getReducedCount(offset);

        return representativeCount;
    }

}