package org.broadinstitute.sting.playground.utils;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.BasicGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.GenotypeGenerator;
import org.broadinstitute.sting.utils.genotype.calls.GenotypeCall;
import org.broadinstitute.sting.utils.genotype.calls.SSGGenotypeCall;
import org.broadinstitute.sting.utils.genotype.confidence.BayesianConfidenceScore;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class GenotypeLikelihoods implements GenotypeGenerator {
    // precalculate these for performance (pow/log10 is expensive!)

    /**
     * SOLID data uses Q0 bases to represent reference-fixed bases -- they shouldn't be counted
     * during GL calculations.  If this field is true, Q0 bases will be removed in add().
     */
    private boolean filterQ0Bases = true;

    private static final double[] oneMinusData = new double[Byte.MAX_VALUE];
    private static final double[] oneHalfMinusDataArachne = new double[Byte.MAX_VALUE];
    private static final double[] oneHalfMinusData3Base = new double[Byte.MAX_VALUE];
    private final boolean keepQ0Bases;
    private static final double log10Of1_3 = log10(1.0 / 3);
    private static final double log10Of2_3 = log10(2.0 / 3);

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            oneMinusData[qual] = log10(1.0 - pow(10, (qual / -10.0)));
        }
    }

    private static double getOneMinusQual(final byte qual) {
        return oneMinusData[qual];
    }


    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            double e = pow(10, (qual / -10.0));
            oneHalfMinusDataArachne[qual] = log10(0.5 - e / 2.0);
            oneHalfMinusData3Base[qual] = log10(0.5 - e / 2.0 + e / 6.0);
            //System.out.printf("qual=%d, e=%f, oneHalfMinusDataArachne=%f, oneHalfMinusData3Base=%f%n", qual, e, oneHalfMinusDataArachne[qual], oneHalfMinusData3Base[qual]);
        }
    }

    private double getOneHalfMinusQual(final byte qual) {
        return oneHalfMinusData[qual];
    }

    public double[] likelihoods;
    public static String[] genotypes = new String[10];
    static {
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
    public int coverage;

    // The genotype priors;
    private double priorHomRef;
    private double priorHet;
    private double priorHomVar;
    private double[] oneHalfMinusData;
    private boolean threeBaseErrors = false;

    // Store the 2nd-best base priors for on-genotype primary bases
    private HashMap<String, Double> onNextBestBasePriors = new HashMap<String, Double>();

    // Store the 2nd-best base priors for off-genotype primary bases
    private HashMap<String, Double> offNextBestBasePriors = new HashMap<String, Double>();

    public GenotypeLikelihoods() {
        double[] p2ndon = {0.000, 0.302, 0.366, 0.142, 0.000, 0.548, 0.370, 0.000, 0.319, 0.000};
        double[] p2ndoff = {0.480, 0.769, 0.744, 0.538, 0.575, 0.727, 0.768, 0.589, 0.762, 0.505};
        keepQ0Bases = true;
        initialize(false, 1.0 - 1e-3, 1e-3, 1e-5, p2ndon, p2ndoff);
    }

    public GenotypeLikelihoods(boolean threeBaseErrors , double priorHomRef, double priorHet, double priorHomVar) {
        double[] p2ndon = {0.000, 0.302, 0.366, 0.142, 0.000, 0.548, 0.370, 0.000, 0.319, 0.000};
        double[] p2ndoff = {0.480, 0.769, 0.744, 0.538, 0.575, 0.727, 0.768, 0.589, 0.762, 0.505};
        keepQ0Bases = true;
        initialize(threeBaseErrors, priorHomRef, priorHet, priorHomVar, p2ndon, p2ndoff);
    }

    public GenotypeLikelihoods(boolean threeBaseErrors , double priorHomRef, double priorHet, double priorHomVar, double[] p2ndon, double[] p2ndoff, boolean keepQ0Bases) {
        this.keepQ0Bases = keepQ0Bases;
        initialize(threeBaseErrors, priorHomRef, priorHet, priorHomVar, p2ndon, p2ndoff);
    }

    /**
     * Are we ignoring Q0 bases during calculations?
     * @return
     */
    public boolean isFilteringQ0Bases() {
        return filterQ0Bases;
    }

    /**
     * Enable / disable filtering of Q0 bases.  Enabled by default
     *
     * @param filterQ0Bases
     */
    public void filterQ0Bases(boolean filterQ0Bases) {
        this.filterQ0Bases = filterQ0Bases;
    }

    private void initialize(boolean threeBaseErrors , double priorHomRef, double priorHet, double priorHomVar, double[] p2ndon, double[] p2ndoff) {
        this.threeBaseErrors  = threeBaseErrors ;
        this.oneHalfMinusData = threeBaseErrors ? oneHalfMinusData3Base : oneHalfMinusDataArachne;

        this.priorHomRef = priorHomRef;
        this.priorHet = priorHet;
        this.priorHomVar = priorHomVar;

        likelihoods = new double[10];

		coverage = 0;

		for (int i = 0; i < likelihoods.length; i++) { likelihoods[i] = Math.log10(0.1); }

        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
            onNextBestBasePriors.put(genotypes[genotypeIndex], p2ndon[genotypeIndex]);
            offNextBestBasePriors.put(genotypes[genotypeIndex], p2ndoff[genotypeIndex]);
        }
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

    public double[] getOnGenotypeSecondaryPriors() {
        double[] p2ndon = new double[10];

        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
            p2ndon[genotypeIndex] = onNextBestBasePriors.get(genotypes[genotypeIndex]);
        }

        return p2ndon;
    }

    public void setOnGenotypeSecondaryPriors(double[] p2ndon) {
        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
            onNextBestBasePriors.put(genotypes[genotypeIndex], p2ndon[genotypeIndex]);
        }
    }

    public double[] getOffGenotypeSecondaryPriors() {
        double[] p2ndoff = new double[10];

        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
            p2ndoff[genotypeIndex] = offNextBestBasePriors.get(genotypes[genotypeIndex]);
        }

        return p2ndoff;
    }

    public void setOffGenotypeSecondaryPriors(double[] p2ndoff) {
        for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
            offNextBestBasePriors.put(genotypes[genotypeIndex], p2ndoff[genotypeIndex]);
        }
    }

    public int add(char ref, char read, byte qual)
	{
        if (qual <= 0) {
            if ( isFilteringQ0Bases() ) {
                return 0;
            } else {
                qual = 1;
            }
        }

		if (coverage == 0)
		{
			for (int i = 0; i < likelihoods.length; i++)
			{
				likelihoods[i] = 0;
			}
		}

		for (int i = 0; i < genotypes.length; i++)
		{
			double likelihood = calculateAlleleLikelihood(ref, read, genotypes[i], qual);
            //if ( originalQual == 0 ) System.out.printf("Likelihood is %f for %c %c %d %s%n", likelihood, ref, read, qual, genotypes[i]);
            likelihoods[i] += likelihood;
			coverage += 1;
        }

        return 1;
    }

    private double calculateAlleleLikelihood(char ref, char read, String genotype, byte qual) {
        if (qual == 0) {
            // zero quals are wrong
            throw new RuntimeException(String.format("Unexpected Q0 base discovered in calculateAlleleLikelihood: %c %c %d", ref, read, qual));
        }

        char h1 = genotype.charAt(0);
        char h2 = genotype.charAt(1);

        double p_base;

        if ((h1 == h2) && (h1 == read)) {
            // hom
            p_base = getOneMinusQual(qual);
        } else if ( (h1 != h2) && ((h1 == read) || (h2 == read)) ) {
            // het
            p_base = getOneHalfMinusQual(qual);
        } else if ( this.threeBaseErrors ) {
            // error
            //System.out.printf("%s %b %f %f%n", genotype, h1 != h2, log10Of2_3, log10Of1_3 );
            p_base = qual / -10.0 + ( h1 != h2 ? log10Of1_3 : log10Of1_3 );
        } else {
            // error
            p_base = qual / -10.0;
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

    public void ApplyPrior(char ref, double[] allele_likelihoods) {
        int k = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = i; j < 4; j++) {
                if (i == j) {
                    this.likelihoods[k] += Math.log10(allele_likelihoods[i]) + Math.log10(allele_likelihoods[j]);
                } else {
                    this.likelihoods[k] += Math.log10(allele_likelihoods[i]) + Math.log10(allele_likelihoods[j]) + Math.log10(2);
                }
                k++;
            }
        }
        this.sort();
    }

    public void ApplyPrior(char ref) {
        for (int i = 0; i < genotypes.length; i++) {
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
        }
        this.sort();
    }

    public void applySecondBaseDistributionPrior(String primaryBases, String secondaryBases) {
        for (int genotypeIndex = 0; genotypeIndex < genotypes.length; genotypeIndex++) {
            char firstAllele = genotypes[genotypeIndex].charAt(0);
            char secondAllele = genotypes[genotypeIndex].charAt(1);

            int offIsGenotypic = 0;
            int offTotal = 0;

            int onIsGenotypic = 0;
            int onTotal = 0;

            for (int pileupIndex = 0; pileupIndex < primaryBases.length(); pileupIndex++) {
                char primaryBase = primaryBases.charAt(pileupIndex);

                if (secondaryBases != null) {
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
            }

            double offPrior = MathUtils.binomialProbability(offIsGenotypic, offTotal, offNextBestBasePriors.get(genotypes[genotypeIndex]));
            double onPrior = MathUtils.binomialProbability(onIsGenotypic, onTotal, onNextBestBasePriors.get(genotypes[genotypeIndex]));

            double logOffPrior = MathUtils.compareDoubles(offPrior, 0.0, 1e-10) == 0 ? Math.log10(Double.MIN_VALUE) : Math.log10(offPrior);
            double logOnPrior = MathUtils.compareDoubles(onPrior, 0.0, 1e-10) == 0 ? Math.log10(Double.MIN_VALUE) : Math.log10(onPrior);

            likelihoods[genotypeIndex] += logOffPrior + logOnPrior;
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

    public double RefPosterior(char ref) {
        this.LodVsRef(ref);
        return this.ref_likelihood;
    }

    private IndelLikelihood indel_likelihood;

    public void addIndelLikelihood(IndelLikelihood indel_likelihood) {
        this.indel_likelihood = indel_likelihood;
    }

    public IndelLikelihood getIndelLikelihood() {
        return this.indel_likelihood;
    }

    /**
     * given all the data associated with a locus, make a genotypeLocus object containing the likelihoods and posterior probs
     *
     * @param tracker contains the reference meta data for this location, which may contain relevent information like dpSNP or hapmap information
     * @param ref     the reference base
     * @param pileup  a pileup of the reads, containing the reads and their offsets
     *
     * @return a GenotypeLocus, containing each of the genotypes and their associated likelihood and posterior prob values
     */
    @Override
    public GenotypeCall callGenotypes(RefMetaDataTracker tracker, char ref, ReadBackedPileup pileup) {
        //filterQ0Bases(!keepQ0Bases); // Set the filtering / keeping flag


        // for calculating the rms of the mapping qualities
        double squared = 0.0;
        for (int i = 0; i < pileup.getReads().size(); i++) {
            SAMRecord read = pileup.getReads().get(i);
            squared += read.getMappingQuality() * read.getMappingQuality();
            int offset = pileup.getOffsets().get(i);
            char base = read.getReadString().charAt(offset);
            byte qual = read.getBaseQualities()[offset];
            add(ref, base, qual);            
        }
        // save off the likelihoods
        if (likelihoods == null || likelihoods.length == 0) return null;

        double lklihoods[] = new double[likelihoods.length];

        System.arraycopy(likelihoods, 0, lklihoods, 0, likelihoods.length);
        

        ApplyPrior(ref);

        applySecondBaseDistributionPrior(pileup.getBases(), pileup.getSecondaryBasePileup());

        // lets setup the locus
        List<Genotype> lst = new ArrayList<Genotype>();
        for (int x = 0; x < this.likelihoods.length; x++) {
                lst.add(new BasicGenotype(pileup.getLocation(),this.genotypes[x],new BayesianConfidenceScore(this.likelihoods[x])));
        }
        return new SSGGenotypeCall(ref,2,pileup.getLocation(),lst,likelihoods,pileup);
    }
    //TODO: add this to the above code
    boolean discovery = false;
    public void setDiscovery(boolean isInDiscoveryMode) {
       discovery = isInDiscoveryMode; 
    }
}
