package org.broadinstitute.sting.utils.genotype.vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class VCFGenotypeRecord
 *         <p/>
 *         The genotype record in VCF store a considerable amount of information,
 *         so they were broken off into their own class
 */
public class VCFGenotypeRecord {
    // the symbol for a empty genotype
    public static final String EMPTY_GENOTYPE = ".";

    // what kind of phasing this genotype has
    public enum PHASE {
        UNPHASED, PHASED, PHASED_SWITCH_PROB, UNKNOWN
    }

    // our phasing
    private PHASE mPhaseType;

    // our bases(s)
    private final List<VCFGenotypeEncoding> mGenotypeAlleles = new ArrayList<VCFGenotypeEncoding>();

    // our mapping of the format mFields to values
    private final Map<String, String> mFields = new HashMap<String, String>();

    // our sample name
    private final String mSampleName;

    /**
     * Create a VCF genotype record
     *
     * @param sampleName
     * @param genotypes
     * @param phasing
     * @param otherFlags
     */
    public VCFGenotypeRecord(String sampleName, List<VCFGenotypeEncoding> genotypes, PHASE phasing, Map<String, String> otherFlags) {
        this.mSampleName = sampleName;
        if (genotypes != null) this.mGenotypeAlleles.addAll(genotypes);
        this.mPhaseType = phasing;
        if (otherFlags != null) this.mFields.putAll(otherFlags);        
    }


    /**
     * determine the phase of the genotype
     *
     * @param phase the string that contains the phase character
     */
    static PHASE determinePhase(String phase) {
        // find the phasing information
        if (phase.equals("/"))
            return PHASE.UNPHASED;
        else if (phase.equals("|"))
            return PHASE.PHASED;
        else if (phase.equals("\\"))
            return PHASE.PHASED_SWITCH_PROB;
        else
            throw new IllegalArgumentException("Unknown genotype phasing parameter");
    }

    /** getter methods */

    public PHASE getPhaseType() {
        return mPhaseType;
    }

    public String getSampleName() {
        return mSampleName;
    }

    public List<VCFGenotypeEncoding> getAlleles() {
        return mGenotypeAlleles;
    }

    public Map<String, String> getFields() {
        return mFields;
    }

    public String toGenotypeString(List<VCFGenotypeEncoding> altAlleles) {
        String str = "";
        boolean first = true;
        for (VCFGenotypeEncoding allele : mGenotypeAlleles) {
            str += String.valueOf((altAlleles.contains(allele)) ? altAlleles.indexOf(allele) + 1 : 0);
            if (first) {
                switch (mPhaseType) {
                    case UNPHASED:
                        str += "/";
                        break;
                    case PHASED:
                        str += "|";
                        break;
                    case PHASED_SWITCH_PROB:
                        str += "\\";
                        break;
                    case UNKNOWN:
                        throw new UnsupportedOperationException("Unknown phase type");
                }
                first = false;
            }
        }        
        return str;

    }


    public boolean equals(Object other) {
        if (other instanceof VCFGenotypeRecord) {
            if (((VCFGenotypeRecord) other).mPhaseType != this.mPhaseType) return false;
            if (!((VCFGenotypeRecord) other).mGenotypeAlleles.equals(this.mGenotypeAlleles)) return false;
            if (!((VCFGenotypeRecord) other).mFields.equals(mFields)) return false;
            if (!((VCFGenotypeRecord) other).mSampleName.equals(this.mSampleName)) return false;
            return true;
        }
        return false;
    }
}
