package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.*;

/**
 * @author ebanks
 *         <p/>
 *         Class VCFGenotypeLocusData
 *         <p/>
 *         represents the meta data for a genotype object.
 */
public class VCFGenotypeLocusData implements GenotypeLocusData, ConfidenceBacked, SLODBacked, AlleleFrequencyBacked {

    // the discovery lod score
    private double mConfidence = 0.0;

    // the strand score lod
    private double mSLOD = 0.0;

    // the allele frequency
    private double mAlleleFrequency = 0.0;

    // the location
    private GenomeLoc mLoc;

    // the ref base
    private char mRefBase;

    /**
     * create a basic genotype meta data pbject, given the following fields
     *
     * @param ref       the reference base
     * @param loc       the locus
     */
    public VCFGenotypeLocusData(char ref, GenomeLoc loc) {
        mRefBase = ref;
        mLoc = loc;
    }

    /**
      * get the reference base.
      * @return a character, representing the reference base
      */
    public char getReference() {
        return mRefBase;
    }

    /**
     * get the genotype's location
     *
     * @return a GenomeLoc representing the location
     */
    public GenomeLoc getLocation() {
        return mLoc;
    }

    /**
     * get the confidence
     *
     * @return the confidence
     */
    public double getConfidence() {
        return mConfidence;
    }

    /**
     *
     * @param   confidence    the confidence for this genotype
     */
    public void setConfidence(double confidence) {
        mConfidence = confidence;
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
     *
     * @param   slod    the strand lod for this genotype
     */
    public void setSLOD(double slod) {
        mSLOD = slod;
    }

    /**
     * get the allele frequency
     *
     * @return the allele frequency
     */
    public double getAlleleFrequency() {
        return mAlleleFrequency;
    }

    /**
     *
     * @param   frequency    the allele frequency for this genotype
     */
    public void setAlleleFrequency(double frequency) {
        mAlleleFrequency = frequency;
    }
}