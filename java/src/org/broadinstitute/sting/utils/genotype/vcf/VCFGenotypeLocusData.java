package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.HashMap;
import java.util.Map;

/**
 * @author ebanks
 *         <p/>
 *         Class VCFGenotypeLocusData
 *         <p/>
 *         represents the meta data for a genotype object.
 */
public class VCFGenotypeLocusData implements GenotypeLocusData, ConfidenceBacked, SLODBacked, IDBacked, AlleleFrequencyBacked, ArbitraryFieldsBacked {

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

    // the id
    private String mID;

    // the various info field values
    private Map<String, String> mInfoFields;

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

    /**
     * @return returns the dbsnp id for this genotype
     */
    public String getID() {
        return mID;
    }

    public void setID(String id) {
        mID = id;
    }

    /**
     * @return returns te arbitrary info fields
     */
    public Map<String, String> getFields() {
        return mInfoFields;
    }

    public void setFields(Map<String, String> fields) {
        mInfoFields = new HashMap<String, String>(fields);
    }    
}