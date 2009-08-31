package org.broadinstitute.sting.playground.gatk.walkers.variants;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;

import java.util.List;

public class IVFBinomialStrand implements IndependentVariantFeature {
    private double strandBalance = 0.5;

    private double[] likelihoods;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            strandBalance = Double.valueOf(arguments);
        }
    }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        likelihoods = new double[10];
        
        ReadBackedPileup pileup = new ReadBackedPileup(context.getReferenceContext().getBase(), context.getAlignmentContext());
        List<SAMRecord> reads = context.getAlignmentContext().getReads();
        String bases = pileup.getBases();

        for (int genotypeIndex = 0; genotypeIndex < Genotype.values().length; genotypeIndex++) {
            Genotype genotype = Genotype.values()[genotypeIndex];

            char[] alleles = { genotype.toString().charAt(0), genotype.toString().charAt(1) };
            int[] strandCounts = { 0, 0 };

            for (int pileupIndex = 0; pileupIndex < bases.length(); pileupIndex++) {
                for (int alleleIndex = 0; alleleIndex < alleles.length; alleleIndex++) {
                    if (bases.charAt(pileupIndex) == alleles[alleleIndex]) {
                        if (reads.get(pileupIndex).getReadNegativeStrandFlag()) {
                            strandCounts[1]++;
                        } else {
                            strandCounts[0]++;
                        }
                    }
                }
            }

            likelihoods[genotypeIndex] = Math.log10(MathUtils.binomialProbability(strandCounts[0], strandCounts[0] + strandCounts[1], strandBalance));
        }
    }

    public double[] getLikelihoods() {
        return likelihoods;
    }

    public String getStudyHeader() {
        return "";
    }

    public String getStudyInfo() {
        return "";
    }
}
