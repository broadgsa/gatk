package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.BaseUtils;

import static java.lang.Math.log10;
import java.util.TreeMap;
import java.util.EnumMap;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

public class ThreeStateErrorGenotypeLikelihoods extends GenotypeLikelihoods {
    //
    // forwarding constructors -- don't do anything at all
    //
    public ThreeStateErrorGenotypeLikelihoods() { super(); }
    public ThreeStateErrorGenotypeLikelihoods(DiploidGenotypePriors priors) { super(priors); }

    /**
     * Cloning of the object
     * @return
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    /**
     * Simple log10(3) cached value
     */
    protected static final double log103 = log10(3.0);

    protected double log10PofTrueBaseGivenMiscall(char observedBase, char chromBase, SAMRecord read, int offset) {
        return -log103; // equivalent to e / 3 model
    }
}
