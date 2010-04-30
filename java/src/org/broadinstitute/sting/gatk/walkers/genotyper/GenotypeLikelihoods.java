package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

/**
 * Stable, error checking version of the Bayesian genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(G | D) = P(G) * P(D | G)
 *
 * where
 *
 * P(D | G) = sum_i log10 P(bi | G)
 *
 * and
 *
 * P(bi | G) = 1 - P(error | q1) if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for homozygous genotypes and for heterozygous genotypes:
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for each of the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 *
 * The priors contain the relative probabilities of each genotype, and must be provided at object creation.
 * From then on, you can call any of the add() routines to update the likelihoods and posteriors in the above
 * model.
 */
public class GenotypeLikelihoods implements Cloneable {
    protected final static int FIXED_PLOIDY = 2;
    protected final static int MAX_PLOIDY = FIXED_PLOIDY + 1;

    protected boolean enableCacheFlag = true;
    protected boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] log10Likelihoods = null;
    protected double[] log10Posteriors = null;

    protected DiploidGenotypePriors priors = null;

    protected FourBaseProbabilities fourBaseLikelihoods = null;

    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     *
     * @param m base model
     */
    public GenotypeLikelihoods(BaseMismatchModel m) {
        this.priors = new DiploidGenotypePriors();
        initialize(m, null);
    }

    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     *
     * @param m base model
     * @param pl default platform
     */
    public GenotypeLikelihoods(BaseMismatchModel m, EmpiricalSubstitutionProbabilities.SequencerPlatform pl) {
        this.priors = new DiploidGenotypePriors();
        initialize(m, pl);
    }

    /**
     * Create a new GenotypeLikelhoods object with given priors for each diploid genotype
     *
     * @param m base model
     * @param priors priors
     */
    public GenotypeLikelihoods(BaseMismatchModel m, DiploidGenotypePriors priors) {
        this.priors = priors;
        initialize(m, null);
    }

    /**
     * Create a new GenotypeLikelhoods object with given priors for each diploid genotype
     *
     * @param m base model
     * @param priors priors
     * @param pl default platform
     */
    public GenotypeLikelihoods(BaseMismatchModel m, DiploidGenotypePriors priors, EmpiricalSubstitutionProbabilities.SequencerPlatform pl) {
        this.priors = priors;
        initialize(m, pl);
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        GenotypeLikelihoods c = (GenotypeLikelihoods)super.clone();
        c.priors = priors;
        c.log10Likelihoods = log10Likelihoods.clone();
        c.log10Posteriors = log10Posteriors.clone();
        c.fourBaseLikelihoods = (FourBaseProbabilities)fourBaseLikelihoods.clone();
        return c;
    }

    protected void initialize(BaseMismatchModel m, EmpiricalSubstitutionProbabilities.SequencerPlatform pl) {
        fourBaseLikelihoods = FourBaseProbabilitiesFactory.makeFourBaseLikelihoods(m, pl);
        setToZero();
    }

    protected void setToZero() {
        log10Likelihoods = zeros.clone();                 // likelihoods are all zeros
        log10Posteriors = priors.getPriors().clone();     // posteriors are all the priors
    }

    public void setVerbose(boolean v) {
        VERBOSE = v;
        fourBaseLikelihoods.setVerbose(v);
    }

    public boolean isVerbose() {
        return VERBOSE;
    }

    public int getMinQScoreToInclude() {
        return fourBaseLikelihoods.getMinQScoreToInclude();
    }

    public void setMinQScoreToInclude(int minQScoreToInclude) {
        fourBaseLikelihoods.setMinQScoreToInclude(minQScoreToInclude);
    }

    /**
     * Returns an array of log10 likelihoods for each genotype, indexed by DiploidGenotype.ordinal values()
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        return log10Likelihoods;
    }

    /**
     * Returns the likelihood associated with DiploidGenotype g
     * @param g genotype
     * @return log10 likelihood as a double
     */
    public double getLikelihood(DiploidGenotype g) {
        return getLikelihoods()[g.ordinal()];
    }

    /**
     * Returns an array of posteriors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return raw log10 (not-normalized posteriors) as a double array
     */
    public double[] getPosteriors() {
        return log10Posteriors;
    }

    /**
     * Returns the posterior associated with DiploidGenotype g
     * @param g genotpe
     * @return raw log10 (not-normalized posteror) as a double
     */
    public double getPosterior(DiploidGenotype g) {
        return getPosteriors()[g.ordinal()];
    }


    /**
     * Returns an array of posteriors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return normalized posterors as a double array
     */
    public double[] getNormalizedPosteriors() {
        double[] normalized = new double[log10Posteriors.length];
        double sum = 0.0;

        // for precision purposes, we need to add (or really subtract, since everything is negative)
        // the largest posterior value from all entries so that numbers don't get too small
        double maxValue = log10Posteriors[0];
        for (int i = 1; i < log10Posteriors.length; i++) {
            if ( maxValue < log10Posteriors[i] )
                maxValue = log10Posteriors[i];
        }

        // collect the posteriors
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double posterior = Math.pow(10, getPosterior(g) - maxValue);
            normalized[g.ordinal()] = posterior;
            sum += posterior;
        }

        // normalize
        for (int i = 0; i < normalized.length; i++)
            normalized[i] /= sum;

        return normalized;
    }



    public DiploidGenotypePriors getPriorObject() {
        return priors;
    }

    /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors.getPriors();
    }

    /**
     * Sets the priors
     * @param priors priors
     */
    public void setPriors(DiploidGenotypePriors priors) {
        this.priors = priors;
        log10Posteriors = zeros.clone();
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Posteriors[i] = priors.getPriors()[i] + log10Likelihoods[i];
        }
    }

    /**
     * Returns the prior associated with DiploidGenotype g
     * @param g genotype
     * @return log10 prior as a double
     */
    public double getPrior(DiploidGenotype g) {
        return getPriors()[g.ordinal()];
    }

    /**
     * Simple function to overload to control the caching of genotype likelihood outputs.
     * By default the function trues true -- do enable caching.  If you are experimenting with an
     * complex calcluation of P(B | G) and caching might not work correctly for you, overload this
     * function and return false, so the super() object won't try to cache your GL calculations.
     *
     * @return true if caching should be enabled, false otherwise
     */
    public boolean cacheIsEnabled() {
        return enableCacheFlag;
    }

    public void setEnableCacheFlag(boolean enable) {
        enableCacheFlag = enable;
    }

    public int add(ReadBackedPileup pileup) {
        return add(pileup, false, false);
    }

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup                    read pileup
     * @param ignoreBadBases            should we ignore bad bases?
     * @param capBaseQualsAtMappingQual should we cap a base's quality by its read's mapping quality?
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual) {
        int n = 0;

        for ( PileupElement p : pileup ) {
            // ignore deletions
            if ( p.isDeletion() )
                continue;

            char base = (char)p.getBase();
            if ( ! ignoreBadBases || ! badBase(base) ) {
                byte qual = capBaseQualsAtMappingQual ? (byte)Math.min((int)p.getQual(), p.getMappingQual()) : p.getQual();
                n += add(base, qual, p.getRead(), p.getOffset());
            }
        }
        
        return n;
    }

    public int add(char observedBase, byte qualityScore, SAMRecord read, int offset) {

        // Handle caching if requested.  Just look up the cached result if its available, or compute and store it
        GenotypeLikelihoods gl;
        if ( cacheIsEnabled() ) {
            if ( ! inCache( observedBase, qualityScore, FIXED_PLOIDY, read) ) {
                 gl = calculateCachedGenotypeLikelihoods(observedBase, qualityScore, FIXED_PLOIDY, read, offset);
            } else {
                gl = getCachedGenotypeLikelihoods(observedBase, qualityScore, FIXED_PLOIDY, read);
            }
        } else {
            gl = calculateGenotypeLikelihoods(observedBase, qualityScore, read, offset);
        }

        // for bad bases, there are no likelihoods
        if ( gl == null )
            return 0;

        double[] likelihoods = gl.getLikelihoods();

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihood = likelihoods[g.ordinal()];
            
            if ( VERBOSE ) {
                boolean fwdStrand = ! read.getReadNegativeStrandFlag();
                System.out.printf("  L(%c | G=%s, Q=%d, S=%s) = %f / %f%n",
                        observedBase, g, qualityScore, fwdStrand ? "+" : "-", pow(10,likelihood) * 100, likelihood);
            }

            log10Likelihoods[g.ordinal()] += likelihood;
            log10Posteriors[g.ordinal()] += likelihood;
        }

        return 1;
    }

    static GenotypeLikelihoods[][][][][][] CACHE = new GenotypeLikelihoods[BaseMismatchModel.values().length][EmpiricalSubstitutionProbabilities.SequencerPlatform.values().length][BaseUtils.BASES.length][QualityUtils.MAX_QUAL_SCORE+1][MAX_PLOIDY][2];

    protected boolean inCache( char observedBase, byte qualityScore, int ploidy, SAMRecord read) {
        return getCache(CACHE, observedBase, qualityScore, ploidy, read) != null;
    }

    protected GenotypeLikelihoods getCachedGenotypeLikelihoods( char observedBase, byte qualityScore, int ploidy, SAMRecord read) {
        GenotypeLikelihoods gl = getCache(CACHE, observedBase, qualityScore, ploidy, read);
        if ( gl == null )
            throw new RuntimeException(String.format("BUG: trying to fetch an unset cached genotype likelihood at base=%c, qual=%d, ploidy=%d, read=%s",
                    observedBase, qualityScore, ploidy, read));
        return gl;
    }

    protected GenotypeLikelihoods calculateCachedGenotypeLikelihoods(char observedBase, byte qualityScore, int ploidy, SAMRecord read, int offset) {
        GenotypeLikelihoods gl = calculateGenotypeLikelihoods(observedBase, qualityScore, read, offset);
        setCache(CACHE, observedBase, qualityScore, ploidy, read, gl);
        return gl;
    }

    protected void setCache( GenotypeLikelihoods[][][][][][] cache,
                             char observedBase, byte qualityScore, int ploidy,
                             SAMRecord read, GenotypeLikelihoods val ) {
        int m = FourBaseProbabilitiesFactory.getBaseMismatchModel(fourBaseLikelihoods).ordinal();
        int a = fourBaseLikelihoods.getReadSequencerPlatformIndex(read);
        int i = BaseUtils.simpleBaseToBaseIndex(observedBase);
        int j = qualityScore;
        int k = ploidy;
        int x = strandIndex(! read.getReadNegativeStrandFlag() );

        cache[m][a][i][j][k][x] = val;
    }

    protected GenotypeLikelihoods getCache( GenotypeLikelihoods[][][][][][] cache,
                                            char observedBase, byte qualityScore, int ploidy, SAMRecord read) {
        int m = FourBaseProbabilitiesFactory.getBaseMismatchModel(fourBaseLikelihoods).ordinal();
        int a = fourBaseLikelihoods.getReadSequencerPlatformIndex(read);
        int i = BaseUtils.simpleBaseToBaseIndex(observedBase);
        int j = qualityScore;
        int k = ploidy;
        int x = strandIndex(! read.getReadNegativeStrandFlag() );
        return cache[m][a][i][j][k][x];
    }

    protected GenotypeLikelihoods calculateGenotypeLikelihoods(char observedBase, byte qualityScore, SAMRecord read, int offset) {
        FourBaseProbabilities fbl = fourBaseLikelihoods.computeLog10Likelihoods(observedBase, qualityScore, read, offset);
        if ( fbl == null )
            return null;

        double[] fbLikelihoods = fbl.getLog10Likelihoods();
        try {

            GenotypeLikelihoods gl = (GenotypeLikelihoods)this.clone();
            gl.setToZero();

            // we need to adjust for ploidy.  We take the raw p(obs | chrom) / ploidy, which is -log10(ploidy) in log space
            double ploidyAdjustment = log10(FIXED_PLOIDY);

            for ( DiploidGenotype g : DiploidGenotype.values() ) {

                // todo assumes ploidy is 2 -- should be generalized.  Obviously the below code can be turned into a loop
                double p_base = 0.0;
                p_base += pow(10, fbLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base1)] - ploidyAdjustment);
                p_base += pow(10, fbLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base2)] - ploidyAdjustment);
                double likelihood = log10(p_base);

                gl.log10Likelihoods[g.ordinal()] += likelihood;
                gl.log10Posteriors[g.ordinal()] += likelihood;
            }

            if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
            }

            return gl;

         } catch ( CloneNotSupportedException e ) {
             throw new RuntimeException(e);
         }
    }


    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // helper routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    public static int strandIndex(boolean fwdStrand) {
        return fwdStrand ? 0 : 1;
    }

    /**
     * Returns true when the observedBase is considered bad and shouldn't be processed by this object.  A base
     * is considered bad if:
     *
     *   Criterion 1: observed base isn't a A,C,T,G or lower case equivalent
     *
     * @param observedBase observed base
     * @return true if the base is a bad base
     */
    protected boolean badBase(char observedBase) {
        return BaseUtils.simpleBaseToBaseIndex(observedBase) == -1;
    }

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        double sum = 0;
        StringBuilder s = new StringBuilder();
        for (DiploidGenotype g : DiploidGenotype.values()) {
            s.append(String.format("%s %.10f ", g, log10Likelihoods[g.ordinal()]));
			sum += Math.pow(10,log10Likelihoods[g.ordinal()]);
        }
		s.append(String.format(" %f", sum));
        return s.toString();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Validation routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(boolean throwException) {
        try {
            priors.validate(throwException);

            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                String bad = null;

                int i = g.ordinal();
                if ( ! MathUtils.wellFormedDouble(log10Likelihoods[i]) || ! MathUtils.isNegativeOrZero(log10Likelihoods[i]) ) {
                    bad = String.format("Likelihood %f is badly formed", log10Likelihoods[i]);
                } else if ( ! MathUtils.wellFormedDouble(log10Posteriors[i]) || ! MathUtils.isNegativeOrZero(log10Posteriors[i]) ) {
                    bad = String.format("Posterior %f is badly formed", log10Posteriors[i]);
                }

                if ( bad != null ) {
                    throw new IllegalStateException(String.format("At %s: %s", g.toString(), bad));
                }
            }
        } catch ( IllegalStateException e ) {
            if ( throwException )
                throw new RuntimeException(e);
            else
                return false;
        }

        return true;
    }

    //
    // Constant static data
    //
    protected final static double[] zeros = new double[DiploidGenotype.values().length];

    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            zeros[g.ordinal()] = 0.0;
        }
    }
}
