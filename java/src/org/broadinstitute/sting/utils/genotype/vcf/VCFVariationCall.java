package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

/**
 * @author ebanks
 *         <p/>
 *         Class VCFVariationCall
 *         <p/>
 *         represents a VCF Variation
 */
public class VCFVariationCall implements VariationCall, VariantBackedByGenotype, ConfidenceBacked, SLODBacked, IDBacked, ArbitraryFieldsBacked {

    // the discovery lod score
    private double mConfidence = 0.0;

    // the strand score lod
    private Double mSLOD = null;

    // the allele frequency field
    private double mAlleleFrequency = -1.0;

    // the location
    private GenomeLoc mLoc;

    // the ref base and alt bases
    private char mRefBase;
    private List<String> mAltBases = new ArrayList<String>();

    // the variant type
    private VARIANT_TYPE mType = VARIANT_TYPE.SNP;

    // the id
    private String mID;

    // the genotypes
    private List<Genotype> mGenotypes = null;

    // the various info field values
    private Map<String, String> mInfoFields;

    /**
     * create a VCF Variation object, given the following fields
     *
     * @param ref       the reference base
     * @param loc       the locus
     * @param type      the variant type
     */
    public VCFVariationCall(char ref, GenomeLoc loc, VARIANT_TYPE type) {
        mRefBase = ref;
        mLoc = loc;
        mType = type;
    }

    /**
      * get the reference base.
      * @return a character, representing the reference base
      */
    public String getReference() {
        return String.valueOf(mRefBase);
    }

    /**
     * get the genotype's location
     *
     * @return a GenomeLoc representing the location
     */
    public GenomeLoc getLocation() {
        return mLoc;
    }

    public boolean isBiallelic() {
        return mAltBases.size() == 1;
    }

    public boolean isSNP() {
        return mType == VARIANT_TYPE.SNP;
    }

    public boolean isInsertion() {
        return mType == VARIANT_TYPE.INSERTION;
    }

    public boolean isIndel() {
        return mType == VARIANT_TYPE.INSERTION || mType == VARIANT_TYPE.DELETION;
    }

    public boolean isDeletion() {
        return mType == VARIANT_TYPE.DELETION;
    }

    public boolean isReference() {
        return mType == VARIANT_TYPE.REFERENCE;
    }

    public VARIANT_TYPE getType() {
        return mType;
    }

    public boolean hasNonRefAlleleFrequency() {
        return mAlleleFrequency >= 0.0;
    }

    public double getNonRefAlleleFrequency() {
        return mAlleleFrequency;
    }

    public double getNegLog10PError() {
        return mConfidence / 10.0;
    }

    public List<String> getAlternateAlleleList() {
        return mAltBases;
    }

    public void addAlternateAllele(String alt) {
        mAltBases.add(alt);
    }

    public List<String> getAlleleList() {
        LinkedList<String> alleles = new LinkedList<String>(mAltBases);
        alleles.addFirst(getReference());
        return alleles;
    }

    public char getAlternativeBaseForSNP() {
        if ( !isSNP() )
            throw new IllegalStateException("This variant is not a SNP");
        if ( mAltBases.size() == 0 )
            throw new IllegalStateException("No alternate alleles have been set");
        return mAltBases.get(0).charAt(0);
    }

    public char getReferenceForSNP() {
        if ( !isSNP() )
            throw new IllegalStateException("This variant is not a SNP");
        return mRefBase;
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
    public Double getSLOD() {
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
     *
     * @param   frequency    the allele frequency for this genotype
     */
    public void setNonRefAlleleFrequency(double frequency) {
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
     *
     * @param   calls    the GenotypeCalls for this variation
     */
    public void setGenotypeCalls(List<Genotype> calls) {
        mGenotypes = calls;
    }

    /**
     * @return a specific genotype that represents the called genotype
     */
    public Genotype getCalledGenotype() {
        if ( mGenotypes == null || mGenotypes.size() != 1 )
            throw new IllegalStateException("There is not one and only one Genotype associated with this Variation");
        return mGenotypes.get(0);
    }

    /**
     * @return a list of all the genotypes
     */
    public List<Genotype> getGenotypes() {
        return mGenotypes;
    }

    /**
     * do we have the specified genotype?  not all backedByGenotypes
     * have all the genotype data.
     *
     * @param x the genotype
     *
     * @return true if available, false otherwise
     */
    public boolean hasGenotype(DiploidGenotype x) {
        if ( mGenotypes == null )
            return false;

        for ( Genotype g : mGenotypes ) {
            if ( DiploidGenotype.valueOf(g.getBases()).equals(x) )
                return true;
        }

        return false;
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