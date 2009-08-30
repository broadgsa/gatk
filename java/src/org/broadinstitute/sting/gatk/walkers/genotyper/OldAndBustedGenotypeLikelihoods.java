package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.BasicGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.confidence.BayesianConfidenceScore;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class OldAndBustedGenotypeLikelihoods {
    protected static final double[] oneMinusData = new double[Byte.MAX_VALUE];
    protected static final double[] oneHalfMinusDataArachne = new double[Byte.MAX_VALUE];
    protected static final double[] oneHalfMinusData3Base = new double[Byte.MAX_VALUE];
    //protected static final double[] oneHalfMinusData = new double[Byte.MAX_VALUE];
    protected static final double log10Of1_3 = log10(1.0 / 3.0);
    private boolean filterQ0Bases = true;

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            double e = pow(10, (qual / -10.0));
            oneMinusData[qual] = log10(1.0 - e);
            oneHalfMinusDataArachne[qual] = log10(0.5 - e / 2.0);
            oneHalfMinusData3Base[qual] = log10(0.5 - e / 2.0 + e / 6.0);
        }
    }

    private static double getOneMinusQual(final byte qual) {
        return oneMinusData[qual];
    }

    private double getOneHalfMinusQual(final byte qual) {
        return oneHalfMinusData[qual];
    }

    //public double[] likelihoods;
    public int coverage;

    // The genotype priors;
    private double priorHomRef;
    private double priorHet;
    private double priorHomVar;
    private double[] oneHalfMinusData;
    public double[] likelihoods;

    public boolean isThreeStateErrors() {
        return threeStateErrors;
    }

    public void setThreeStateErrors(boolean threeStateErrors) {
        this.threeStateErrors = threeStateErrors;
        this.oneHalfMinusData = threeStateErrors ? oneHalfMinusData3Base : oneHalfMinusDataArachne;
    }

    private boolean threeStateErrors = false;
    private boolean discoveryMode = false;

    public static double[] computePriors(double h) {
        double[] pdbls = new double[3];
        pdbls[0] = 1.0 - (3.0 * h / 2.0);
        pdbls[1] = h;
        pdbls[2] = h / 2.0;
        return pdbls;
    }

    public final static String[] genotypes = new String[10];
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

    /**
     * set the mode to discovery
     * @param isInDiscoveryMode
     */
    public void setDiscovery(boolean isInDiscoveryMode) {
       discoveryMode = isInDiscoveryMode;
    }
    // Store the 2nd-best base priors for on-genotype primary bases
    private HashMap<String, Double> onNextBestBasePriors = new HashMap<String, Double>();

    // Store the 2nd-best base priors for off-genotype primary bases
    private HashMap<String, Double> offNextBestBasePriors = new HashMap<String, Double>();

    public OldAndBustedGenotypeLikelihoods(double heterozygosity) {
        double[] vals = computePriors(heterozygosity);
        initialize(vals[0], vals[1], vals[2]);
    }

    public OldAndBustedGenotypeLikelihoods(double priorHomRef, double priorHet, double priorHomVar) {
        initialize(priorHomRef, priorHet, priorHomVar);
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

    private void initialize(double priorHomRef, double priorHet, double priorHomVar) {
        this.oneHalfMinusData = threeStateErrors ? oneHalfMinusData3Base : oneHalfMinusDataArachne;

        this.priorHomRef = priorHomRef;
        this.priorHet = priorHet;
        this.priorHomVar = priorHomVar;

        likelihoods = new double[10];

		coverage = 0;

		for (int i = 0; i < likelihoods.length; i++) { likelihoods[i] = Math.log10(0.1); }
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
        //System.out.printf("%s %s%n", OldAndBustedGenotypeLikelihoods.class, this.toString('A'));

		for (int i = 0; i < genotypes.length; i++)
		{
			double likelihood = calculateAlleleLikelihood(ref, read, genotypes[i], qual);
            System.out.printf("Likelihood is %f for %c %c %d %s%n", likelihood, ref, read, qual, genotypes[i]);
            likelihoods[i] += likelihood;
			coverage += 1;
        }

        //System.out.printf("%s %s%n", OldAndBustedGenotypeLikelihoods.class, this.toString('A'));

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
        } else if ( threeStateErrors ) {
            // error
            //System.out.printf("%s %c %c %b %f %f%n", genotype, h1, h2, h1 != h2, log10Of2_3, log10Of1_3 );
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
            s = s + sorted_genotypes[i] + ":" + String.format("%.10f", sorted_likelihoods[i]);
			sum += Math.pow(10,sorted_likelihoods[i]);
        }
		s = s + String.format(" %f", sum);
        return s;
    }

    public void applyPrior(char ref, double[] allele_likelihoods) {
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

    public void applyPrior(char ref) {
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

    /**
     * given all the data associated with a locus, make a genotypeLocus object containing the likelihoods and posterior probs
     *
     * @param tracker contains the reference meta data for this location, which may contain relevent information like dpSNP or hapmap information
     * @param ref     the reference base
     * @param pileup  a pileup of the reads, containing the reads and their offsets
     *
     * @return a GenotypeLocus, containing each of the genotypes and their associated likelihood and posterior prob values
     */
    public SSGGenotypeCall callGenotypes(RefMetaDataTracker tracker, char ref, ReadBackedPileup pileup) {
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

        // Apply the two calculations
        applyPrior(ref);

        // lets setup the locus
        List<Genotype> lst = new ArrayList<Genotype>();
        for (int x = 0; x < this.likelihoods.length; x++) {
                lst.add(new BasicGenotype(pileup.getLocation(),this.genotypes[x],new BayesianConfidenceScore(this.likelihoods[x])));
        }
        return new SSGGenotypeCall(discoveryMode,ref,2,lst,likelihoods,pileup);
    }
}