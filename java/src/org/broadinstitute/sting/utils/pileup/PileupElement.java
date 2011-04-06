package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.utils.*;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class PileupElement {
    public static final byte DELETION_BASE = BaseUtils.D;
    public static final byte DELETION_QUAL = (byte) 16;
    public static final byte A_FOLLOWED_BY_INSERTION_BASE = (byte) 87;
    public static final byte C_FOLLOWED_BY_INSERTION_BASE = (byte) 88;
    public static final byte T_FOLLOWED_BY_INSERTION_BASE = (byte) 89;
    public static final byte G_FOLLOWED_BY_INSERTION_BASE = (byte) 90;

    protected SAMRecord read;
    protected int offset;

    public PileupElement( SAMRecord read, int offset ) {
        this.read = read;
        this.offset = offset;
    }

    public boolean isDeletion() {
        return offset == -1;
    }

    public SAMRecord getRead() { return read; }
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

    public int getMappingQual() { return read.getMappingQuality(); }

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
}