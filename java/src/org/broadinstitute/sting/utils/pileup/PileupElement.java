package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class PileupElement {
    public static final byte DELETION_BASE = BaseUtils.D;
    public static final byte DELETION_QUAL = 0;

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

    public byte getSecondBase() {
        return getSecondBase(offset);
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
        return isDeletion() ? DELETION_BASE : BaseUtils.simpleBaseToBaseIndex((char)read.getReadBases()[offset]);
    }

    protected byte getSecondBase(final int offset) {
        return isDeletion() ? DELETION_BASE : BaseUtils.getSecondBase(read, offset);
    }

    protected byte getQual(final int offset) {
        return isDeletion() ? DELETION_QUAL : read.getBaseQualities()[offset];
    }
}