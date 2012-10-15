package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

@Deprecated
public class ConstrainedDiploidExactAFCalc extends DiploidExactAFCalc {
    protected ConstrainedDiploidExactAFCalc(int nSamples, int maxAltAlleles, int maxAltAllelesForIndels, final int ploidy) {
        super(nSamples, maxAltAlleles, maxAltAllelesForIndels, ploidy);
    }

    protected StateTracker makeMaxLikelihood(final VariantContext vc, final AFCalcResultTracker resultTracker) {
        final int[] maxACsToConsider = computeMaxACs(vc);
        resultTracker.setAClimits(maxACsToConsider);
        return new StateTracker(maxACsToConsider);
    }

    /**
     * Computes the maximum ACs we need to consider for each alt allele
     *
     * Walks over the genotypes in VC, and computes for each alt allele the maximum
     * AC we need to consider in that alt allele dimension.  Does the calculation
     * based on the PLs in each genotype g, choosing to update the max AC for the
     * alt alleles corresponding to that PL.  Only takes the first lowest PL,
     * if there are multiple genotype configurations with the same PL value.  It
     * takes values in the order of the alt alleles.
     *
     * @param vc the variant context we will compute max alt alleles for
     * @return a vector of max alt alleles, indexed by alt allele, so result[0] is the AC of the
     *          first alt allele.
     */
    @Ensures("result != null")
    protected final int[] computeMaxACs(final VariantContext vc) {
        final int[] maxACs = new int[vc.getNAlleles()-1];

        for ( final Genotype g : vc.getGenotypes() )
            updateMaxACs(g, maxACs);

        return maxACs;
    }

    /**
     * Update the maximum achievable allele counts in maxAC according to the PLs in g
     *
     * Selects the maximum genotype configuration from the PLs in g, and updates
     * the maxAC for this configure.  For example, if the lowest PL is for 0/1, updates
     * the maxAC for the alt allele 1 by 1.  If it's 1/1, update is 2.  Works for
     * many number of alt alleles (determined by length of maxACs).
     *
     * If the max PL occurs at 0/0, updates nothing
     * Note that this function greedily takes the first min PL, so that if 0/1 and 1/1 have
     * the same PL value, then updates the first one.
     *
     * Also, only will update 1 alt allele, so if 0/1 and 0/2 both have the same PL,
     * then only first one (1) will be updated
     *
     * @param g the genotype to update
     * @param maxACs the max allele count vector for alt alleles (starting at 0 => first alt allele)
     */
    @Requires({
            "g != null",
            "maxACs != null",
            "goodMaxACs(maxACs)"})
    private void updateMaxACs(final Genotype g, final int[] maxACs) {
        final int[] PLs = g.getLikelihoods().getAsPLs();

        int minPLi = 0;
        int minPL = PLs[0];

        for ( int i = 0; i < PLs.length; i++ ) {
            if ( PLs[i] < minPL ) {
                minPL = PLs[i];
                minPLi = i;
            }
        }

        final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair pair = GenotypeLikelihoods.getAllelePair(minPLi);
        updateMaxACs(maxACs, pair.alleleIndex1);
        updateMaxACs(maxACs, pair.alleleIndex2);
    }

    /**
     * Simple helper.  Update max alt alleles maxACs according to the allele index (where 0 == ref)
     *
     * If alleleI == 0 => doesn't update anything
     * else maxACs[alleleI - 1]++
     *
     * @param maxACs array of max alt allele ACs
     * @param alleleI the index (relative to 0) to update a count of 1 in max alt alleles.
     */
    @Requires({
            "alleleI >= 0",
            "(alleleI - 1) < maxACs.length",
            "goodMaxACs(maxACs)"})
    private void updateMaxACs(final int[] maxACs, final int alleleI) {
        if ( alleleI > 0 )
            maxACs[alleleI-1]++;
    }

    private static boolean goodMaxACs(final int[] maxACs) {
        return MathUtils.sum(maxACs) >= 0;
    }
}
