package org.broadinstitute.sting.utils.pileup;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class ExtendedPileupElement extends PileupElement {
    private int pileupOffset = 0;
    private ReadBackedPileup pileup = null;

    public ExtendedPileupElement( SAMRecord read, int readOffset, int pileupOffset, ReadBackedPileup pileup ) {
        super(read, readOffset);
        this.pileupOffset = pileupOffset;
        this.pileup = pileup;
    }

    public int getPileupOffset() {
        return pileupOffset;
    }

    public ReadBackedPileup getPileup() {
        return pileup;
    }
}