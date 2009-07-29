package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class GenotypeBucket
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public class GenotypeLocusImpl implements GenotypeLocus {

    private final List<Genotype> mGenotypes = new ArrayList<Genotype>();
    private GenomeLoc mLocation = null;
    private int mReadDepth = -1;
    private double mRMSMappingQual = -1;

    public GenotypeLocusImpl(GenomeLoc location, int readDepth, double rmsMappingQual) {
        this.mLocation = location;
        mReadDepth = readDepth;
        mRMSMappingQual = rmsMappingQual;
    }

    /**
     * Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     * is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    @Override
    public GenomeLoc getLocation() {
        return mLocation;
    }

    /**
     * get the ploidy at this locus
     *
     * @return an integer representing the genotype ploidy at this location
     */
    @Override
    public int getPloidy() {
        return 2;
    }

    /**
     * get the genotypes, sorted in asscending order by their likelihoods (the best
     * to the worst likelihoods)
     *
     * @return a list of the likelihoods
     */
    @Override
    public List<Genotype> getGenotypes() {
        Collections.sort(this.mGenotypes, new LikelihoodComparator());
        return this.mGenotypes;
    }

    

    /**
     * get the genotypes and their posteriors
     *
     * @return a list of the poseriors
     */
    @Override
    public List<Genotype> getPosteriors() {
        Collections.sort(this.mGenotypes, new PosteriorComparator());
        return this.mGenotypes;
    }

    /**
     * get the genotypes sorted lexigraphically
     *
     * @return a list of the genotypes sorted lexi
     */
    @Override
    public List<Genotype> getLexigraphicallySortedGenotypes() {
        Collections.sort(this.mGenotypes, new LexigraphicalComparator());
        return this.mGenotypes;
    }

    /**
     * get the read depth at this position
     *
     * @return the read depth, -1 if it is unknown
     */
    @Override
    public int getReadDepth() {
        return mReadDepth;
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
        this.mGenotypes.add(genotype);
    }

    /**
     * get the root mean square (RMS) of the mapping qualities
     *
     * @return the RMS, or a value < 0 if it's not available
     */
    @Override
    public double getRMSMappingQuals() {
        return mRMSMappingQual;
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
