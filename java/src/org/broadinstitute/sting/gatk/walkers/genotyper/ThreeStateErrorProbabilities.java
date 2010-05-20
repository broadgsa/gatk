package org.broadinstitute.sting.gatk.walkers.genotyper;

import static java.lang.Math.log10;

import net.sf.samtools.SAMRecord;

public class ThreeStateErrorProbabilities extends FourBaseProbabilities {
    //
    // forwarding constructors -- don't do anything at all
    //
    public ThreeStateErrorProbabilities() { super(); }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    /**
     * Simple log10(3) cached value
     */
    protected static final double log103 = log10(3.0);

    /**
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param read         SAM read
     * @param offset       offset on read
     * @return log10 likelihood
     */
    protected double log10PofTrueBaseGivenMiscall(byte observedBase, byte chromBase, SAMRecord read, int offset) {
        return -log103; // equivalent to e / 3 model
    }
}