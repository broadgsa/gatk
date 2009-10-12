package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import net.sf.samtools.SAMRecord;

/**
 * This implements the old style CRD calculation of the chance that a base being a true chromBase given
 * an miscalled base, in which the p is e, grabbing all of the probability.  It shouldn't be used
 */
public class OneStateErrorGenotypeLikelihoods extends GenotypeLikelihoods {
    //
    // forwarding constructors -- don't do anything at all
    //
    public OneStateErrorGenotypeLikelihoods() { super(); }
    public OneStateErrorGenotypeLikelihoods(DiploidGenotypePriors priors) { super(priors); }

    /**
     * Cloning of the object
     * @return
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    /**
     *
     * @param observedBase
     * @param chromBase
     * @param read
     * @param offset
     * @return
     */
    protected double log10PofTrueBaseGivenMiscall(char observedBase, char chromBase, SAMRecord read, int offset) {
        return 0; // equivalent to e model
    }

    static GenotypeLikelihoods[][][][] ONE_STATE_CACHE = new GenotypeLikelihoods[BaseUtils.BASES.length][QualityUtils.MAX_QUAL_SCORE][MAX_PLOIDY][2];
    protected GenotypeLikelihoods getSetCache( char observedBase, byte qualityScore, int ploidy,
                                             SAMRecord read, int offset, GenotypeLikelihoods val ) {
        return simpleGetSetCache(ONE_STATE_CACHE, observedBase, qualityScore, ploidy, read, offset, val);
    }
}
