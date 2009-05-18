package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.MathUtils;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import java.util.HashMap;

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

    // Store the 2nd-best base priors for on-genotype primary bases
    HashMap<String, Double> onNextBestBasePriors = new HashMap<String, Double>();

    // Store the 2nd-best base priors for off-genotype primary bases
    HashMap<String, Double> offNextBestBasePriors = new HashMap<String, Double>();

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

        onNextBestBasePriors.put("AA", 0.000);
        onNextBestBasePriors.put("AC", 0.302);
        onNextBestBasePriors.put("AG", 0.366);
        onNextBestBasePriors.put("AT", 0.142);
        onNextBestBasePriors.put("CC", 0.000);
        onNextBestBasePriors.put("CG", 0.548);
        onNextBestBasePriors.put("CT", 0.370);
        onNextBestBasePriors.put("GG", 0.000);
        onNextBestBasePriors.put("GT", 0.319);
        onNextBestBasePriors.put("TT", 0.000);

        offNextBestBasePriors.put("AA", 0.480);
        offNextBestBasePriors.put("AC", 0.769);
        offNextBestBasePriors.put("AG", 0.744);
        offNextBestBasePriors.put("AT", 0.538);
        offNextBestBasePriors.put("CC", 0.575);
        offNextBestBasePriors.put("CG", 0.727);
        offNextBestBasePriors.put("CT", 0.768);
        offNextBestBasePriors.put("GG", 0.589);
        offNextBestBasePriors.put("GT", 0.762);
        offNextBestBasePriors.put("TT", 0.505);
    }

    public void add(char ref, char read, byte qual) {
        for (int i = 0; i < genotypes.length; i++) {
            likelihoods[i] += calculateAlleleLikelihood(ref, read, genotypes[i], qual);
        }
    }

    private double calculateAlleleLikelihood(char ref, char read, String genotype, byte qual) 
	{
		if (qual == 0) { qual = 1; } // zero quals are wrong.

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

    public void ApplyPrior(char ref, char alt, double p_alt) 
    {
        for (int i = 0; i < genotypes.length; i++) 
		{
            if ((p_alt == -1) || (p_alt <= 1e-6))
            {
	            if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) 
				{
	                // hom-ref
	                likelihoods[i] += Math.log10(1.0 - 1e-3);
	            } 
				else if ((genotypes[i].charAt(0) != ref) && (genotypes[i].charAt(1) != ref)) 
				{
	                // hom-nonref
	                likelihoods[i] += Math.log10(1e-5);
	            } 
				else 
				{
	                // het
	                likelihoods[i] += Math.log10(1e-3);
	            }
                if (Double.isInfinite(likelihoods[i])) { likelihoods[i] = -1000; }
            }
            else
            {
	            if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref)) 
				{
	                // hom-ref
	                likelihoods[i] += 2.0 * Math.log10(1.0 - p_alt);
	            } 
				else if ((genotypes[i].charAt(0) == alt) && (genotypes[i].charAt(1) == alt)) 
				{
	                // hom-nonref
	                likelihoods[i] += 2.0 * Math.log10(p_alt);
	            } 
				else if (((genotypes[i].charAt(0) == alt) && (genotypes[i].charAt(1) == ref)) ||
						 ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == alt)))
				{
	                // het
	                likelihoods[i] += Math.log10((1.0-p_alt) * p_alt * 2.0);
	            }
				else
				{
					// something else (noise!)
					likelihoods[i] += Math.log10(1e-5);
				}

                if (Double.isInfinite(likelihoods[i])) { likelihoods[i] = -1000; }
            }
        }
        this.sort();
    }

    public void applyFourBaseDistributionPrior(String primaryBases, String secondaryBases) {
        for (int genotypeIndex = 0; genotypeIndex < genotypes.length; genotypeIndex++) {
            char firstAllele = genotypes[genotypeIndex].charAt(0);
            char secondAllele = genotypes[genotypeIndex].charAt(1);

            int offIsGenotypic = 0;
            int offTotal = 0;

            int onIsGenotypic = 0;
            int onTotal = 0;

            for (int pileupIndex = 0; pileupIndex < primaryBases.length(); pileupIndex++) {
                char primaryBase = primaryBases.charAt(pileupIndex);
                char secondaryBase = secondaryBases.charAt(pileupIndex);

                if (primaryBase != firstAllele && primaryBase != secondAllele) {
                    if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
                        offIsGenotypic++;
                    }
                    offTotal++;
                } else {
                    if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
                        onIsGenotypic++;
                    }
                    onTotal++;
                }
            }

            double offPrior = MathUtils.binomialProbability(offIsGenotypic, offTotal, offNextBestBasePriors.get(genotypes[genotypeIndex]));
            double onPrior = MathUtils.binomialProbability(onIsGenotypic, onTotal, onNextBestBasePriors.get(genotypes[genotypeIndex]));

            likelihoods[genotypeIndex] += Math.log10(offPrior) + Math.log10(onPrior);

            //System.out.println(genotypes[genotypeIndex] + " " + offNextBestBasePriors.get(genotypes[genotypeIndex]) + " " + offIsGenotypic + " " + offTotal + " " + (((double) offIsGenotypic)/((double) offTotal)) + " " + offPrior);
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

    public AlleleFrequencyEstimate toAlleleFrequencyEstimate(GenomeLoc location, char ref, int depth, String bases, double[] posteriors, String sample_name) {
        this.sort();
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
            qhat = 1.0;
            qstar = 1.0;
            alt = sorted_genotypes[0].charAt(0);
        } else {
            // het
            qhat = 0.5;
            qstar = 0.5;

            if (sorted_genotypes[0].charAt(0) != ref) {
                alt = sorted_genotypes[0].charAt(0);
            }
            if (sorted_genotypes[0].charAt(1) != ref) {
                alt = sorted_genotypes[0].charAt(1);
            }
        }

		this.LodVsRef(ref); //HACK
		//System.out.printf("DBG: %f %f\n", sorted_likelihoods[0], ref_likelihood); 

        return new AlleleFrequencyEstimate(location, ref, alt, 2, qhat, qstar, this.LodVsRef(ref), this.LodVsNextBest(), sorted_likelihoods[0], ref_likelihood, depth, bases, (double[][]) null, this.likelihoods, sample_name);
    }

	public static class IndelCall
	{
		public String   type;
		public String[] alleles;
		public double   p;
		public double   lod;

		public IndelCall(String type, String[] alleles, double p, double lod)
		{
			this.type = type;
			this.alleles = alleles;
			this.p = p;
			this.lod = lod;
		}

		public String toString()
		{
			return String.format("%s/%s %f %f", alleles[0], alleles[1], p, lod); 
		}
	}

	public static IndelCall callIndel(String[] indels)
	{
		HashMap<String,Integer> indel_allele_counts = new HashMap<String,Integer>();
		for (int i = 0; i < indels.length; i++)
		{
			if (! indel_allele_counts.containsKey(indels[i])) 
			{
				indel_allele_counts.put(indels[i], 1);
			}
			else
			{
				indel_allele_counts.put(indels[i], indel_allele_counts.get(indels[i])+1);
			}
		}

		Object[] keys = indel_allele_counts.keySet().toArray();
		String[] alleles = new String[keys.length];
		int[] counts = new int[keys.length];
		double likelihoods[] = new double[keys.length];
		int null_count = 0;
		String max_allele = null;
		int max_count = -1;
		if ((keys.length > 0) && (! ((keys.length == 1) && (((String)keys[0]).equals("null")))))
		{
			for (int i = 0; i < keys.length; i++)
			{
				Integer count = (Integer)indel_allele_counts.get(keys[i]);
				alleles[i] = (String)keys[i];
				counts[i] = count;
				if (alleles[i].equals("null")) { null_count = counts[i]; }
				else if (counts[i] > max_count) { max_count = counts[i]; max_allele = alleles[i]; }
				//System.out.printf("%s[%d] ", keys[i], count);
			}
			//System.out.printf("\n");

			double eps = 1e-3;
			double pRef = null_count*Math.log10(1.0 - eps)   + max_count*Math.log10(eps) + Math.log10(0.999);
			double pHet = null_count*Math.log10(0.5 - eps/2) + max_count*Math.log10(0.5-eps/2) + Math.log10(1e-3);
			double pHom = null_count*Math.log10(eps)         + max_count*Math.log10(1.0 - eps) + Math.log10(1e-5);

			double lodRef = pRef - Math.max(pHet, pHom);
			double lodHet = pHet - pRef;
			double lodHom = pHom - pRef;

			//System.out.printf("%s/%s %f %f\n", "null", "null", pRef, lodRef);
			//System.out.printf("%s/%s %f %f\n", max_allele, "null", pHet, lodHet);
			//System.out.printf("%s/%s %f %f\n", max_allele, max_allele, pHom, lodHom);
			//System.out.printf("\n");

			if (lodRef > 0)
			{
				// reference call
				String[] genotype = new String[2];
				genotype[0] = "null";
				genotype[1] = "null";
				return new IndelCall("ref", genotype, pRef, lodRef);
			}
			else if (lodHet > lodHom)
			{
				// het call
				String[] genotype = new String[2];
				genotype[0] = "null";
				genotype[1] = max_allele;
				return new IndelCall("het", genotype, pHet, lodHet);
			}
			else
			{
				// hom call
				String[] genotype = new String[2];
				genotype[0] = max_allele;
				genotype[1] = max_allele;
				return new IndelCall("hom", genotype, pHom, lodHom);
			}
		}
		return null;
	}

}
