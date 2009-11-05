package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;

/**
 * This implements the old style CRD calculation of the chance that a base being a true chromBase given
 * an miscalled base, in which the p is e, grabbing all of the probability.  It shouldn't be used
 */
public class OneStateErrorProbabilities extends FourBaseProbabilities {
    //
    // forwarding constructors -- don't do anything at all
    //
    public OneStateErrorProbabilities() { super(); }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    /**
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param read         SAM read
     * @param offset       offset on read
     * @return log10 likelihood
     */
    protected double log10PofTrueBaseGivenMiscall(char observedBase, char chromBase, SAMRecord read, int offset) {
        return 0; // equivalent to e model
    }
}