package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.ArrayList;
import java.util.List;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public class VariantCounter extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    long nBasesCovered = 0;
    int nSNPs = 0;
    int nHets = 0;

    public VariantCounter() {
        super("variant_counts");
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        nSNPs += eval == null ? 0 : 1;

        if ( this.getMaster().evalContainsGenotypes && eval != null && eval.isSNP() && ((VariantBackedByGenotype)eval).getGenotype( DiploidGenotype.valueOf(eval.getAlternateBases())).isHet() )
            nHets++;

        return null;
    }

    /**
     * No need to finalize the data in general
     * @param nSites
     */
    public void finalize(long nSites) {
        nBasesCovered = nSites;
    }

    private double variantRate(int n) {
        return n / (1.0 * Math.max(nBasesCovered, 1));
    }

    private long variantRateInverse(int n) {
        return nBasesCovered / Math.max(n, 1);
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("n bases covered    %d", nBasesCovered));
        s.add(String.format("variants           %d", nSNPs));
        s.add(String.format("variant rate       %.5f confident variants per base", variantRate(nSNPs)));
        s.add(String.format("variant rate       1 / %d confident variants per base [human single sample genome-wide expectation is ~1 / 666]", variantRateInverse(nSNPs)));

        if ( this.getMaster().evalContainsGenotypes ) {
            s.add(String.format("heterozygotes      %d", nHets));
            s.add(String.format("homozygotes        %d", nSNPs - nHets));

            s.add(String.format("heterozygosity     %.5f confident hets per base", variantRate(nHets)));
            s.add(String.format("heterozygosity     1 / %d confident hets per base [human single sample expectation is ~1 / 1000]", variantRateInverse(nHets)));

            s.add(String.format("het to hom ratio   %.2f confident hets per confident homozygote non-refs [human single sample genome-wide expectation is 2:1]",
                    ((double)nHets) / (Math.max(nSNPs - nHets, 1))));
        }
        return s;
    }
}