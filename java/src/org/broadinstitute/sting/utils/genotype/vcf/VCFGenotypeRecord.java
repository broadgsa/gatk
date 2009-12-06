package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class VCFGenotypeRecord
 *         <p/>
 */
public class VCFGenotypeRecord implements Genotype {
    // the symbols for an empty genotype
    public static final String EMPTY_GENOTYPE = "./.";
    public static final String EMPTY_ALLELE = ".";

    public static final int MISSING_DEPTH = -1;

    // what kind of phasing this genotype has
    public enum PHASE {
        UNPHASED, PHASED, PHASED_SWITCH_PROB, UNKNOWN
    }

    // our record
    private VCFRecord mRecord;

    // our phasing
    private PHASE mPhaseType;

    // our bases(s)
    private final List<VCFGenotypeEncoding> mGenotypeAlleles = new ArrayList<VCFGenotypeEncoding>();

    // our mapping of the format mFields to values
    private final Map<String, String> mFields = new HashMap<String, String>();

    // our sample name
    private String mSampleName;

    /**
     * Create a VCF genotype record
     *
     * @param sampleName  sample name
     * @param genotypes   list of genotypes
     * @param phasing     phasing
     * @param otherFlags  other flags
     */
    public VCFGenotypeRecord(String sampleName, List<VCFGenotypeEncoding> genotypes, PHASE phasing, Map<String, String> otherFlags) {
        this.mSampleName = sampleName;
        if (genotypes != null) this.mGenotypeAlleles.addAll(genotypes);
        this.mPhaseType = phasing;
        if (otherFlags != null) this.mFields.putAll(otherFlags);        
    }

    public void setVCFRecord(VCFRecord record) {
        this.mRecord = record;
    }

    public void setSampleName(String name) {
        mSampleName = name;
    }

    /**
     * determine the phase of the genotype
     *
     * @param phase the string that contains the phase character
     *
     * @return the phase
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

    public double getNegLog10PError() {
        return ( mFields.containsKey("GQ") ? Double.valueOf(mFields.get("GQ")) / 10.0 : 0.0);
    }

    public int getReadCount() {
        int depth = MISSING_DEPTH;
        if ( mFields.containsKey("RD") )
            depth = Integer.valueOf(mFields.get("RD"));
        else if ( mFields.containsKey("DP") )
            depth = Integer.valueOf(mFields.get("DP"));
        return depth;
    }

    public GenomeLoc getLocation() {
        return mRecord != null ? mRecord.getLocation() : null;
    }

    public char getReference() {
        return mRecord != null ? mRecord.getReferenceBase() : 'N';
    }

    public Variation toVariation(char ref) {
        return mRecord != null ? mRecord : null;
    }

    public String getBases() {
        String genotype = "";
        for ( VCFGenotypeEncoding encoding : mGenotypeAlleles )
            genotype += encoding.getBases();
        return genotype;
    }

    public boolean isVariant(char ref) {
        for ( VCFGenotypeEncoding encoding : mGenotypeAlleles ) {
            if ( encoding.getType() == VCFGenotypeEncoding.TYPE.UNCALLED )
                continue;
            if ( encoding.getType() != VCFGenotypeEncoding.TYPE.SINGLE_BASE ||
                 encoding.getBases().charAt(0) != ref )
                return true;
        }
        return false;
    }

    public boolean isPointGenotype() {
        return (mRecord != null ? !mRecord.isIndel() : true);
    }

    public boolean isHom() {
        if ( mGenotypeAlleles.size() == 0 )
            return true;

        String bases = mGenotypeAlleles.get(0).getBases();
        for ( int i = 1; i < mGenotypeAlleles.size(); i++ ) {
            if ( !bases.equals(mGenotypeAlleles.get(1).getBases()) )
                return false;
        }
        return true;
    }

    public boolean isHet() {
        return !isHom();
    }

    public int getPloidy() {
        return 2;
    }

    private String toGenotypeString(List<VCFGenotypeEncoding> altAlleles) {
        String str = "";
        boolean first = true;
        for (VCFGenotypeEncoding allele : mGenotypeAlleles) {
            if (allele.getType() == VCFGenotypeEncoding.TYPE.UNCALLED)
                str += VCFGenotypeRecord.EMPTY_ALLELE;
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
        for (String field : mFields.keySet()) {
            if (mFields.get(field).equals("")) continue;
            builder.append(VCFRecord.GENOTYPE_FIELD_SEPERATOR);
            builder.append(mFields.get(field));

        }
        return builder.toString();
    }
}
