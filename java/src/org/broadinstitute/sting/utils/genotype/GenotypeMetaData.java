package org.broadinstitute.sting.utils.genotype;


/**
 * @author ebanks
 *         <p/>
 *         Class GenotypeMetaData
 *         <p/>
 *         represents the meta data for a genotype object.
 */
public class GenotypeMetaData {

    // the discovery lod score
    private double mLOD;

    // the strand score lod
    private double mSLOD;

    // the allele frequency
    private double mAlleleFrequency;

    /**
     * create a basic genotype meta data pbject, given the following fields
     *
     * @param discoveryLOD       the discovery lod
     * @param strandLOD          the strand score lod
     * @param alleleFrequency    the allele frequency
     */
    public GenotypeMetaData(double discoveryLOD, double strandLOD, double alleleFrequency) {
        mLOD = discoveryLOD;
        mSLOD = strandLOD;
        mAlleleFrequency = alleleFrequency;
    }

    /**
     * get the discovery lod
     *
     * @return the discovery lod
     */
    public double getLOD() {
        return mLOD;
    }

    /**
     * get the strand lod
     *
     * @return the strand lod
     */
    public double getSLOD() {
        return mSLOD;
    }

    /**
     * get the allele frequency
     *
     * @return the allele frequency
     */
    public double getAlleleFrequency() {
        return mAlleleFrequency;
    }
}