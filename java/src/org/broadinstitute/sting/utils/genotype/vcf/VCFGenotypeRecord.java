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

    private String toGenotypeString(List<VCFGenotypeEncoding> altAlleles) {
        String str = "";
        boolean first = true;
        for (VCFGenotypeEncoding allele : mGenotypeAlleles) {
            if (allele.getType() == VCFGenotypeEncoding.TYPE.UNCALLED)
                str += VCFGenotypeRecord.EMPTY_GENOTYPE;
            else
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

    public boolean isEmptyGenotype() {
        for ( VCFGenotypeEncoding encoding : mGenotypeAlleles ) {
            if ( encoding.getType() != VCFGenotypeEncoding.TYPE.UNCALLED )
                return false;
        }
        return true;
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

    /**
     * output a string representation of the VCFGenotypeRecord, given the alternate alleles
     *
     * @param altAlleles the alternate alleles, needed for toGenotypeString()
     *
     * @return a string
     */
    public String toStringEncoding(List<VCFGenotypeEncoding> altAlleles) {
        StringBuilder builder = new StringBuilder();
        builder.append(toGenotypeString(altAlleles));
        boolean first = true;
        for (String field : mFields.keySet()) {
            if (mFields.get(field).equals("")) continue;
            builder.append(VCFRecord.GENOTYPE_FIELD_SEPERATOR);
            builder.append(mFields.get(field));

        }
        return builder.toString();
    }
}
