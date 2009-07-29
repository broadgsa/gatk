package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class GenotypeCallImpl
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public class GenotypeCallImpl implements GenotypeCall {

    // our stored genotype locus
    private final GenotypeLocus mLocus;
    private final char mRefBase;
    private final ConfidenceScore mScore;

    /**
     * generate a GenotypeCall object with the specified locus info and reference base
     *
     * @param mLocus   the locus
     * @param mRefBase the reference base to use
     */
    public GenotypeCallImpl(GenotypeLocus mLocus, char mRefBase, ConfidenceScore mScore) {
        if (mLocus.getGenotypes().size() < 1) throw new StingException("Genotype Locus is empty");
        this.mLocus = mLocus;
        this.mRefBase = String.valueOf(mRefBase).toUpperCase().charAt(0);
        this.mScore = mScore;
    }

    /**
     * get the confidence
     *
     * @return a ConfidenceScore representing the LOD score that this genotype was called with
     */
    @Override
    public ConfidenceScore getConfidence() {
        return mScore;
    }

    /**
     * gets the reference base
     *
     * @return the reference base we represent
     */
    @Override
    public char getReferencebase() {
        return mRefBase;
    }

    /**
     * get the best vrs the next best genotype LOD score
     *
     * @return the genotype, and a LOD for best - next
     */
    @Override
    public Pair<Genotype, ConfidenceScore> getBestVrsNext() {
        List<Genotype> genos = this.mLocus.getGenotypes();
        if (mLocus.getGenotypes().size() < 2) throw new StingException("Genotype Locus does not contain two genotypes");
        return new Pair<Genotype, ConfidenceScore>(genos.get(0),
                new ConfidenceScore(Math.abs(genos.get(0).getLikelihood() - genos.get(1).getLikelihood()), ConfidenceScore.SCORE_METHOD.BEST_NEXT));
    }

    /**
     * get the best vrs the reference allele.
     *
     * @return the genotype, and a LOD for best - ref.  The best may be ref, unless you've checked
     *         with is variation
     */
    @Override
    public Pair<Genotype, ConfidenceScore> getBestVrsRef() {
        List<Genotype> genos = this.mLocus.getGenotypes();

        // find the reference allele
        String ref = Utils.dupString(this.mRefBase, mLocus.getPloidy()).toUpperCase();
        Genotype refGenotype = findRefGenotype(ref, genos);
        if (mLocus.getGenotypes().size() < 2) throw new StingException("Genotype Locus does not contain two genotypes");
        return new Pair<Genotype, ConfidenceScore>(genos.get(0),
                new ConfidenceScore(Math.abs(genos.get(0).getLikelihood() - refGenotype.getLikelihood()), ConfidenceScore.SCORE_METHOD.BEST_NEXT));
    }

    /**
     * get the reference genotype object
     *
     * @param ref   the reference as a ploidy count homozygous string
     * @param genos the genotype list
     *
     * @return a genotype for the
     */
    private static Genotype findRefGenotype(String ref, List<Genotype> genos) {
        Genotype refGenotype = null;
        for (Genotype g : genos) {
            if (g.getBases().equals(ref)) refGenotype = g;
        }
        if (refGenotype == null) {
            for (Genotype g : genos) {
                System.err.println(g.getBases());
            }
            throw new StingException("Unable to find the reference genotype + " + ref + " size of genotype list = " + genos.size());
        }
        return refGenotype;
    }

    /**
     * check to see if this call is a variant, i.e. not homozygous reference
     *
     * @return true if it's not hom ref, false otherwise
     */
    @Override
    public boolean isVariation() {
        List<Genotype> genos = this.mLocus.getGenotypes();
        String ref = Utils.dupString(this.mRefBase, mLocus.getPloidy()).toUpperCase();
        return !(genos.get(0).getBases().equals(ref));
    }

    /** return genotype locus, with our data */
    @Override
    public GenotypeLocus toGenotypeLocus() {
        return mLocus;
    }

    /**
     * Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     * is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    @Override
    public GenomeLoc getLocation() {
        return mLocus.getLocation();
    }

    /**
     * get the ploidy at this locus
     *
     * @return an integer representing the genotype ploidy at this location
     */
    @Override
    public int getPloidy() {
        return mLocus.getPloidy();
    }

    /**
     * get the genotypes, sorted in asscending order by their likelihoods (the best
     * to the worst likelihoods)
     *
     * @return a list of the likelihoods
     */
    @Override
    public List<Genotype> getGenotypes() {
        return mLocus.getGenotypes();
    }

    /**
     * get the genotypes and their posteriors
     *
     * @return a list of the poseriors
     */
    @Override
    public List<Genotype> getPosteriors() {
        return mLocus.getPosteriors();
    }

    /**
     * get the genotypes sorted lexigraphically
     *
     * @return a list of the genotypes sorted lexi
     */
    @Override
    public List<Genotype> getLexigraphicallySortedGenotypes() {
        return mLocus.getLexigraphicallySortedGenotypes();
    }

    /**
     * get the read depth at this position
     *
     * @return the read depth, -1 if it is unknown
     */
    @Override
    public int getReadDepth() {
        return mLocus.getReadDepth();
    }

    /**
     * add a genotype to the collection
     *
     * @param genotype
     *
     * @throws InvalidGenotypeException
     */
    @Override
    public void addGenotype(Genotype genotype) throws InvalidGenotypeException {
        mLocus.addGenotype(genotype);
    }

    /**
     * get the root mean square (RMS) of the mapping qualities
     *
     * @return the RMS, or a value < 0 if it's not available
     */
    @Override
    public double getRMSMappingQuals() {
        return mLocus.getRMSMappingQuals();
    }

    /**
     * create a variant, given the reference, and a lod score
     *
     * @param refBase the reference base
     * @param score   the threshold to use to determine if it's a variant or not
     *
     * @return a variant object, or null if no genotypes meet the criteria
     */
    @Override
    public Variant toGenotypeCall(char refBase, ConfidenceScore score) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
