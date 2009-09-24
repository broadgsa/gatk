package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.QualityUtils;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import java.util.HashMap;
import java.lang.Cloneable;

public class ClassicGenotypeLikelihoods implements Cloneable {
    // precalculate these for performance (pow/log10 is expensive!)
    private static final double[] oneMinusData = new double[Byte.MAX_VALUE];

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            oneMinusData[qual] = log10(1.0 - pow(10, (qual / -10.0)));
            //oneMinusData[qual] = log10(1.0 - QualityUtils.qualToProb(qual));
        }
    }

    private static double getOneMinusQual(final byte qual) {
        return oneMinusData[qual];
    }

    private static final double[] oneHalfMinusData = new double[Byte.MAX_VALUE];

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            oneHalfMinusData[qual] = log10(0.5 - pow(10, (qual / -10.0)) / 2.0);
            //oneHalfMinusData[qual] = log10(0.5 - QualityUtils.qualToProb(qual) / 2.0);
        }
    }

    private static double getOneHalfMinusQual(final byte qual) {
        return oneHalfMinusData[qual];
    }

    public double[] likelihoods;
    public String[] genotypes;
	public int coverage;

    // The genotype priors;
    private double priorHomRef;
    private double priorHet;
    private double priorHomVar;
    public String[] sorted_genotypes;
    public double[] sorted_likelihoods;
    double ref_likelihood = Double.NaN;
	private IndelLikelihood indel_likelihood;

    // Store the 2nd-best base priors for on-genotype primary bases
    //private HashMap<String, Double> onNextBestBasePriors = new HashMap<String, Double>();

    // Store the 2nd-best base priors for off-genotype primary bases
    //private HashMap<String, Double> offNextBestBasePriors = new HashMap<String, Double>();

    private static double[] p2ndon = {0.000, 0.302, 0.366, 0.142, 0.000, 0.548, 0.370, 0.000, 0.319, 0.000};
    private static double[] p2ndoff = {0.480, 0.769, 0.744, 0.538, 0.575, 0.727, 0.768, 0.589, 0.762, 0.505};

    public ClassicGenotypeLikelihoods() {
        initialize(1.0 - 1e-3, 1e-3, 1e-5, p2ndon, p2ndoff);
    }

    public ClassicGenotypeLikelihoods(boolean foo) {
    }

    public ClassicGenotypeLikelihoods(double priorHomRef, double priorHet, double priorHomVar) {
        initialize(priorHomRef, priorHet, priorHomVar, p2ndon, p2ndoff);
    }

    public ClassicGenotypeLikelihoods(double priorHomRef, double priorHet, double priorHomVar, double[] p2ndon, double[] p2ndoff) {
        initialize(priorHomRef, priorHet, priorHomVar, p2ndon, p2ndoff);
    }

    public ClassicGenotypeLikelihoods clone() {
        ClassicGenotypeLikelihoods c = new ClassicGenotypeLikelihoods(false);
        c.likelihoods = this.likelihoods.clone();
        c.genotypes = this.genotypes.clone();
	    c.coverage = this.coverage;

        // The genotype priors;
	    c.priorHomRef = this.priorHomRef;
	    c.priorHet = this.priorHet;
	    c.priorHomVar = this.priorHomVar;
        //public String[] sorted_genotypes;
        //public double[] sorted_likelihoods;
        //double ref_likelihood = Double.NaN;
	    //private IndelLikelihood indel_likelihood;
	    return c;
    }

    private void initialize(double priorHomRef, double priorHet, double priorHomVar, double[] p2ndon, double[] p2ndoff) {
        this.priorHomRef = priorHomRef;
        this.priorHet = priorHet;
        this.priorHomVar = priorHomVar;

        likelihoods = new double[10];
        genotypes = new String[10];
		coverage = 0;

		for (int i = 0; i < likelihoods.length; i++) { likelihoods[i] = Math.log10(0.1); }

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

        //for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
        //    onNextBestBasePriors.put(genotypes[genotypeIndex], p2ndon[genotypeIndex]);
        //    offNextBestBasePriors.put(genotypes[genotypeIndex], p2ndoff[genotypeIndex]);
        //}
    }

    public double getHomRefPrior() {
        return priorHomRef;
    }

    public void setHomRefPrior(double priorHomRef) {
        this.priorHomRef = priorHomRef;
    }

    public double getHetPrior() {
        return priorHet;
    }

    public void setHetPrior(double priorHet) {
        this.priorHet = priorHet;
    }

    public double getHomVarPrior() {
        return priorHomVar;
    }

    public void setHomVarPrior(double priorHomVar) {
        this.priorHomVar = priorHomVar;
    }

//    public double[] getOnGenotypeSecondaryPriors() {
//        double[] p2ndon = new double[10];
//
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            p2ndon[genotypeIndex] = onNextBestBasePriors.get(genotypes[genotypeIndex]);
//        }
//
//        return p2ndon;
//    }
//
//    public void setOnGenotypeSecondaryPriors(double[] p2ndon) {
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            onNextBestBasePriors.put(genotypes[genotypeIndex], p2ndon[genotypeIndex]);
//        }
//    }
//
//    public double[] getOffGenotypeSecondaryPriors() {
//        double[] p2ndoff = new double[10];
//
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            p2ndoff[genotypeIndex] = offNextBestBasePriors.get(genotypes[genotypeIndex]);
//        }
//
//        return p2ndoff;
//    }
//
//    public void setOffGenotypeSecondaryPriors(double[] p2ndoff) {
//        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
//            offNextBestBasePriors.put(genotypes[genotypeIndex], p2ndoff[genotypeIndex]);
//        }
//    }

    public void add(char ref, char read, byte qual) 
	{ 
		if (qual <= 0) { qual = 1; }

		if (coverage == 0)
		{
			for (int i = 0; i < likelihoods.length; i++)
			{
				likelihoods[i] = 0;
			}
		}
		double sum = 0.0;
        for (int i = 0; i < genotypes.length; i++) 
		{
			double likelihood = calculateAlleleLikelihood(ref, read, genotypes[i], qual); 
            likelihoods[i] += likelihood;
        }
		coverage += 1;
    }

    public void add(char ref, char read, byte qual, ConfusionMatrix matrix, String platform) 
	{ 
		if (qual <= 0) { qual = 1; }
		if (platform == null) { platform = "ILLUMINA"; }
		if (read == 'N') { return; }

		if (coverage == 0)
		{
			for (int i = 0; i < likelihoods.length; i++)
			{
				likelihoods[i] = 0;
			}
		}
		double sum = 0.0;
        for (int i = 0; i < genotypes.length; i++) 
		{
        	char h1 = genotypes[i].charAt(0);
        	char h2 = genotypes[i].charAt(1);

			double p1 = matrix.lookup(platform, read, h1);
			double p2 = matrix.lookup(platform, read, h2);

			double likelihood = calculateAlleleLikelihood(ref, read, genotypes[i], qual, p1, p2); 

			//System.out.printf("DBG: %c %c %s %d %f %f %f\n", ref, read, genotypes[i], qual, p1, p2, likelihood);

            likelihoods[i] += likelihood;
        }
		coverage += 1;
    }

    private double calculateAlleleLikelihood(char ref, char read, String genotype, byte qual) {
        if (qual == 0) {
            qual = 1;
        } // zero quals are wrong.

        char h1 = genotype.charAt(0);
        char h2 = genotype.charAt(1);

        double p_base;

        if ((h1 == h2) && (h1 == read)) {
            // hom
            p_base = getOneMinusQual(qual);
        } else if ((h1 != h2) && ((h1 == read) || (h2 == read))) {
            // het
            p_base = getOneHalfMinusQual(qual);
        } else {
            // error
            p_base = qual / -10.0;
        }

        return p_base;
    }

	public void TEST()
	{
		double p_A2A = 1.00;
		double p_T2T = 1.00;

		double p_A2T = 0.75;
		double p_T2A = 0.25;

		char ref = 'A';

		System.out.printf("\tA\tT\n");
		System.out.printf("A\t%.02f\t%.02f\n", p_A2A, p_A2T);
		System.out.printf("T\t%.02f\t%.02f\n", p_T2A, p_T2T);
		System.out.printf("\n");

		System.out.printf("P(A,Q20|AA) = %f\n", calculateAlleleLikelihood(ref, 'A', "AA", (byte)20, p_A2A, p_A2A));
		System.out.printf("P(A,Q20|AT) = %f\n", calculateAlleleLikelihood(ref, 'A', "AT", (byte)20, p_A2A, p_A2T));
		System.out.printf("P(A,Q20|TT) = %f\n", calculateAlleleLikelihood(ref, 'A', "TT", (byte)20, p_A2T, p_A2T));

		System.out.printf("P(T,Q20|AA) = %f\n", calculateAlleleLikelihood(ref, 'T', "AA", (byte)20, p_T2A, p_T2A));
		System.out.printf("P(T,Q20|AT) = %f\n", calculateAlleleLikelihood(ref, 'T', "AT", (byte)20, p_T2A, p_T2T));
		System.out.printf("P(T,Q20|TT) = %f\n", calculateAlleleLikelihood(ref, 'T', "TT", (byte)20, p_T2T, p_T2T));

		//System.exit(0);
	}

    private double calculateAlleleLikelihood(char ref, char read, String genotype, byte qual, double p1, double p2) {
        if (qual == 0) {
            qual = 1;
        } // zero quals are wrong.

        char h1 = genotype.charAt(0);
        char h2 = genotype.charAt(1);

		double perr = Math.pow(10.0,qual/-10.0);

        double p_base = 0;

		if (read == h1)
		{
			p_base += (1.0 - perr);
		}
		else
		{
			p_base += perr * p1;
		}

		if (read == h2)
		{
			p_base += (1.0 - perr);
		}
		else
		{
			p_base += perr * p2;
		}

		p_base = Math.log10(p_base/2.0);

        return p_base;
    }

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
		double sum = 0;
        String s = String.format("%s %f %f ", this.BestGenotype(), this.LodVsNextBest(), this.LodVsRef(ref));
        for (int i = 0; i < sorted_genotypes.length; i++) {
            if (i != 0) {
                s = s + " ";
            }
            s = s + sorted_genotypes[i] + ":" + String.format("%.2f", sorted_likelihoods[i]);
			sum += Math.pow(10,sorted_likelihoods[i]);
        }
		s = s + String.format(" %f", sum);
        return s;
    }

    public void ApplyPrior(char ref, double[] allele_likelihoods)
	{
		int k = 0;
		for (int i = 0; i < 4; i++)
		{ 
			for (int j = i; j < 4; j++)
			{
				if (i == j) 
				{
					this.likelihoods[k] += Math.log10(allele_likelihoods[i]) + Math.log10(allele_likelihoods[j]);
				}
				else
				{
					this.likelihoods[k] += Math.log10(allele_likelihoods[i]) + Math.log10(allele_likelihoods[j]) + Math.log10(2);
				}
				k++;
			}
		}
		this.sort();
	}

    public void ApplyPrior(char ref, char alt, double p_alt) {
        for (int i = 0; i < genotypes.length; i++) {
            if ((p_alt == -1) || (p_alt <= 1e-6)) {
                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) {
                    // hom-ref
                    likelihoods[i] += Math.log10(priorHomRef);
                } else if ((genotypes[i].charAt(0) != ref) && (genotypes[i].charAt(1) != ref)) {
                    // hom-nonref
                    likelihoods[i] += Math.log10(priorHomVar);
                } else {
                    // het
                    likelihoods[i] += Math.log10(priorHet);
                }
                if (Double.isInfinite(likelihoods[i])) {
                    likelihoods[i] = -1000;
                }
            } else {
                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) {
                    // hom-ref
                    likelihoods[i] += 2.0 * Math.log10(1.0 - p_alt);
                } else if ((genotypes[i].charAt(0) == alt) && (genotypes[i].charAt(1) == alt)) {
                    // hom-nonref
                    likelihoods[i] += 2.0 * Math.log10(p_alt);
                } else if (((genotypes[i].charAt(0) == alt) && (genotypes[i].charAt(1) == ref)) ||
                        ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == alt))) {
                    // het
                    likelihoods[i] += Math.log10((1.0 - p_alt) * p_alt * 2.0);
                } else {
                    // something else (noise!)
                    likelihoods[i] += Math.log10(1e-5);
                }

                if (Double.isInfinite(likelihoods[i])) {
                    likelihoods[i] = -1000;
                }
            }
        }
        this.sort();
    }

	public void ApplyWeight(double weight)
	{
	    double log10weight = Math.log10(weight);
		for (int i = 0; i < genotypes.length; i++) { likelihoods[i] += log10weight; }
		this.sort();
	}

//    public void applySecondBaseDistributionPrior(String primaryBases, String secondaryBases) {
//        for (int genotypeIndex = 0; genotypeIndex < genotypes.length; genotypeIndex++) {
//            char firstAllele = genotypes[genotypeIndex].charAt(0);
//            char secondAllele = genotypes[genotypeIndex].charAt(1);
//
//            int offIsGenotypic = 0;
//            int offTotal = 0;
//
//            int onIsGenotypic = 0;
//            int onTotal = 0;
//
//            for (int pileupIndex = 0; pileupIndex < primaryBases.length(); pileupIndex++) {
//                char primaryBase = primaryBases.charAt(pileupIndex);
//
//                if (secondaryBases != null) {
//                    char secondaryBase = secondaryBases.charAt(pileupIndex);
//
//                    if (primaryBase != firstAllele && primaryBase != secondAllele) {
//                        if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
//                            offIsGenotypic++;
//                        }
//                        offTotal++;
//                    } else {
//                        if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
//                            onIsGenotypic++;
//                        }
//                        onTotal++;
//                    }
//                }
//            }
//
//            double offPrior = MathUtils.binomialProbability(offIsGenotypic, offTotal, offNextBestBasePriors.get(genotypes[genotypeIndex]));
//            double onPrior = MathUtils.binomialProbability(onIsGenotypic, onTotal, onNextBestBasePriors.get(genotypes[genotypeIndex]));
//
//            likelihoods[genotypeIndex] += Math.log10(offPrior) + Math.log10(onPrior);
//        }
//        this.sort();
//    }

    public double LodVsNextBest() {
        this.sort();
        return sorted_likelihoods[0] - sorted_likelihoods[1];
    }

    public double LodVsRef(char ref) {
        if ((this.BestGenotype().charAt(0) == ref) && (this.BestGenotype().charAt(1) == ref)) {
            ref_likelihood = sorted_likelihoods[0];
            return (-1.0 * this.LodVsNextBest());
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

	public double RefPosterior(char ref)
	{
		this.LodVsRef(ref);
		return this.ref_likelihood;
	}

	public void addIndelLikelihood(IndelLikelihood indel_likelihood) { this.indel_likelihood = indel_likelihood; }
	public IndelLikelihood getIndelLikelihood() { return this.indel_likelihood; }

}
