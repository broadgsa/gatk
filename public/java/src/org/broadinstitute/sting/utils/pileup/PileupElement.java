package org.broadinstitute.sting.utils.pileup;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.BaseUtils;
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

    @Requires({
            "read != null",
            "offset >= -1",
            "offset <= read.getReadLength()"})
    public PileupElement( GATKSAMRecord read, int offset ) {
        this.read = read;
        this.offset = offset;
    }

    public boolean isDeletion() {
        return offset == -1;
    }

    @Ensures("result != null")
    public GATKSAMRecord getRead() { return read; }

    @Ensures("result == offset")
    public int getOffset() { return offset; }

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
        return String.format("%s @ %d = %c Q%d", getRead().getReadName(), getOffset(), (char)getBase(), getQual());
    }

    protected byte getBase(final int offset) {
        return isDeletion() ? DELETION_BASE : read.getReadBases()[offset];
    }

    protected int getBaseIndex(final int offset) {
        return BaseUtils.simpleBaseToBaseIndex(isDeletion() ? DELETION_BASE : read.getReadBases()[offset]);
    }

    protected byte getQual(final int offset) {
        return isDeletion() ? DELETION_QUAL : read.getBaseQualities()[offset];
    }

    @Override
    public int compareTo(final PileupElement pileupElement) {
        if ( offset < pileupElement.offset )
            return -1;
        else if ( offset > pileupElement.offset )
            return 1;
        else if ( read.getAlignmentStart() < pileupElement.read.getAlignmentStart() )
            return -1;
        else if ( read.getAlignmentStart() > pileupElement.read.getAlignmentStart() )
            return 1;
        else
            return 0;
    }

    // --------------------------------------------------------------------------
    //
    // Reduced read accessors
    //
    // --------------------------------------------------------------------------

    public boolean isReducedRead() {
        return ((GATKSAMRecord)read).isReducedRead();
    }

    public int getRepresentativeCount() {
        return isReducedRead() ? ((GATKSAMRecord)read).getReducedCount(offset) : 1;
    }

}