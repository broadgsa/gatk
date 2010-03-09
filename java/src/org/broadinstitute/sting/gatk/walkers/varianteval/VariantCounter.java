package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Genotype;
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
        nSNPs += eval == null || eval.isReference() ? 0 : 1;

        // TODO -- break the het check out to a different module used only for single samples

        if ( this.getMaster().evalContainsGenotypes && eval != null && eval.isBiallelic() && eval.isSNP() && eval instanceof VariantBackedByGenotype ) {
            List<Genotype> genotypes = ((VariantBackedByGenotype)eval).getGenotypes();
            if ( genotypes.size() == 1 && genotypes.get(0).isHet() ) {
                nHets++;
            }
        }

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
        s.add(String.format("variant rate       1 / %d confident variants per base", variantRateInverse(nSNPs)));

        if ( this.getMaster().evalContainsGenotypes ) {
            s.add(String.format("heterozygotes      %d", nHets));
            s.add(String.format("homozygotes        %d", nSNPs - nHets));

            s.add(String.format("heterozygosity     %.5f confident hets per base", variantRate(nHets)));
            s.add(String.format("heterozygosity     1 / %d confident hets per base", variantRateInverse(nHets)));

            s.add(String.format("het to hom ratio   %.2f confident hets per confident homozygote non-refs   [single sample genome-wide expectation is 2:1]",
                    ((double)nHets) / (Math.max(nSNPs - nHets, 1))));
        }
        return s;
    }
}