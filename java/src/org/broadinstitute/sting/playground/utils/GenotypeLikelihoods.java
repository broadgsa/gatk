package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

public class GenotypeLikelihoods {
    // precalculate these for performance (pow/log10 is expensive!)
    private static final double[] oneMinusData = new double[Byte.MAX_VALUE];
    static {
        for(int qual=0; qual < Byte.MAX_VALUE; qual++) {
            oneMinusData[qual] = log10(1.0 - pow(10,(qual/-10.0)));
        }
    }
    private static double getOneMinusQual(final byte qual) {
        return oneMinusData[qual];
    }

    private static final double[] oneHalfMinusData = new double[Byte.MAX_VALUE];
    static {
        for(int qual=0; qual < Byte.MAX_VALUE; qual++) {
            oneHalfMinusData[qual] = log10(0.5-pow(10,(qual/-10.0))/2.0);
        }
    }

    private static double getOneHalfMinusQual(final byte qual) {
        return oneHalfMinusData[qual];
    }



    public double[] likelihoods;
    public String[] genotypes;

    public GenotypeLikelihoods() {
        likelihoods = new double[10];
        genotypes = new String[10];

        genotypes[0] = "AA";
        genotypes[1] = "AC";
        genotypes[2] = "AG";
        genotypes[3] = "AT";
        genotypes[4] = "CC";
        genotypes[5] = "CG";
        genotypes[6] = "CT";
        genotypes[7] = "GG";
        genotypes[8] = "GT";
        genotypes[9] = "TT";
    }

    public void add(char ref, char read, byte qual) {
        for (int i = 0; i < genotypes.length; i++) {
            likelihoods[i] += calculateAlleleLikelihood(ref, read, genotypes[i], qual);
        }
    }

    private double calculateAlleleLikelihood(char ref, char read, String genotype, byte qual) {
        char h1 = genotype.charAt(0);
        char h2 = genotype.charAt(1);

        double p_base;

        if ((h1 == h2) && (h1 == read)) {
            p_base = getOneMinusQual(qual); //Math.log10(1 - p_error);
        } else if ((h1 != h2) && (h1 == read) || (h2 == read)) {
            p_base = getOneHalfMinusQual(qual); // )Math.log10(0.5 - (p_error / 2.0));
        } else {
            // the real math would be
            //     likelihood += log10(pow(10,(qual/-10.0)));
            // but it simplifies to
            p_base = qual/-10.0;
        }

        return p_base;
    }

    public String[] sorted_genotypes;
    public double[] sorted_likelihoods;

    public void sort() {
        Integer[] permutation = Utils.SortPermutation(likelihoods);

        Integer[] reverse_permutation = new Integer[permutation.length];
        for (int i = 0; i < reverse_permutation.length; i++) {
            reverse_permutation[i] = permutation[(permutation.length - 1) - i];
        }

        sorted_genotypes = Utils.PermuteArray(genotypes, reverse_permutation);
        sorted_likelihoods = Utils.PermuteArray(likelihoods, reverse_permutation);
    }

    public String toString(char ref) {
        this.sort();
        String s = String.format("%s %f %f ", this.BestGenotype(), this.LodVsNextBest(), this.LodVsRef(ref));
        for (int i = 0; i < sorted_genotypes.length; i++) {
            if (i != 0) {
                s = s + " ";
            }
            s = s + sorted_genotypes[i] + ":" + String.format("%.2f", sorted_likelihoods[i]);
        }
        return s;
    }

    public void ApplyPrior(char ref, double p_alt) {
        for (int i = 0; i < genotypes.length; i++) {
            if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) {
                // hom-ref
                likelihoods[i] += Math.log10(1.0 - 1e-3);
            } else if ((genotypes[i].charAt(0) != ref) && (genotypes[i].charAt(1) != ref)) {
                // hom-nonref
                likelihoods[i] += Math.log10(1e-5);
            } else {
                // het
                likelihoods[i] += Math.log10(1e-3);
            }
        }
        this.sort();
    }

    public double LodVsNextBest() {
        this.sort();
        return sorted_likelihoods[0] - sorted_likelihoods[1];
    }

    double ref_likelihood = Double.NaN;

    public double LodVsRef(char ref) {
        if ((this.BestGenotype().charAt(0) == ref) && (this.BestGenotype().charAt(1) == ref)) {
            ref_likelihood = sorted_likelihoods[0];
            return this.LodVsNextBest();
        } else {
            for (int i = 0; i < genotypes.length; i++) {
                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) {
                    ref_likelihood = likelihoods[i];
                }
            }
        }
        return sorted_likelihoods[0] - ref_likelihood;
    }

    public String BestGenotype() {
        this.sort();
        return this.sorted_genotypes[0];
    }

    public double BestPosterior() {
        this.sort();
        return this.sorted_likelihoods[0];
    }

    public AlleleFrequencyEstimate toAlleleFrequencyEstimate(GenomeLoc location, char ref, int depth, String bases, double[] posteriors) {
        double qhat = Double.NaN;
        double qstar = Double.NaN;
        char alt = 'N';

        if ((sorted_genotypes[0].charAt(0) == ref) && (sorted_genotypes[0].charAt(1) == ref)) {
            // hom-ref
            qhat = 0.0;
            qstar = 0.0;
            alt = 'N';
        } else if ((sorted_genotypes[0].charAt(0) != ref) && (sorted_genotypes[0].charAt(1) != ref)) {
            // hom-nonref
            likelihoods[0] += Math.log10(1e-5);
            qhat = 1.0;
            qstar = 1.0;
            alt = sorted_genotypes[0].charAt(0);
        } else {
            // het
            likelihoods[0] += Math.log10(1e-3);
            qhat = 0.5;
            qstar = 0.5;

            if (sorted_genotypes[0].charAt(0) != ref) {
                alt = sorted_genotypes[0].charAt(0);
            }
            if (sorted_genotypes[0].charAt(1) != ref) {
                alt = sorted_genotypes[0].charAt(1);
            }
        }

        return new AlleleFrequencyEstimate(location, ref, alt, 2, qhat, qstar, this.LodVsRef(ref), this.LodVsNextBest(), sorted_likelihoods[0], ref_likelihood, depth, bases, (double[][]) null, this.likelihoods);
    }

}
