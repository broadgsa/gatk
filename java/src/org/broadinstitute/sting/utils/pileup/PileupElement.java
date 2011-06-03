package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.utils.*;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import com.google.java.contract.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 */
public class PileupElement {
    public static final byte DELETION_BASE = BaseUtils.D;
    public static final byte DELETION_QUAL = (byte) 16;
    public static final byte A_FOLLOWED_BY_INSERTION_BASE = (byte) 87;
    public static final byte C_FOLLOWED_BY_INSERTION_BASE = (byte) 88;
    public static final byte T_FOLLOWED_BY_INSERTION_BASE = (byte) 89;
    public static final byte G_FOLLOWED_BY_INSERTION_BASE = (byte) 90;

    protected final SAMRecord read;
    protected final int offset;

    @Requires({
            "read != null",
            "offset >= -1",
            "offset <= read.getReadLength()"})
    public PileupElement( SAMRecord read, int offset ) {
        this.read = read;
        this.offset = offset;
    }

    public boolean isDeletion() {
        return offset == -1;
    }

    @Ensures("result != null")
    public SAMRecord getRead() { return read; }

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

    // --------------------------------------------------------------------------
    //
    // Reduced read accessors
    //
    // --------------------------------------------------------------------------

    private Integer getReducedReadQualityTagValue() {
        return (Integer)getRead().getAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG);
    }

    public boolean isReducedRead() {
        return getReducedReadQualityTagValue() != null;
    }

    public int getReducedCount() {
        return (int)getQual();
    }

    public byte getReducedQual() {
        return (byte)(int)getReducedReadQualityTagValue();
    }

}