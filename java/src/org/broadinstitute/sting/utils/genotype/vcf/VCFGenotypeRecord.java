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
    // what kind of phasing this genotype has
    public enum PHASE {
        UNPHASED, PHASED, PHASED_SWITCH_PROB, UNKNOWN
    }

    // our phasing
    private PHASE phaseType;

    // our reference bases(s)
    private final char mReferenceBase;

    // our bases(s)
    private final List<String> mAlleleBases = new ArrayList<String>();

    // our mapping of the format mFields to values
    private final Map<String, String> mFields = new HashMap<String, String>();

    // our sample name
    private final String mSampleName;
    /**
     * create a VCF record
     *
     * @param keyValues     the key values
     * @param Alleles       the alleles, one if we're halpoid, two if we're diploid
     * @param phasing       the phasing of the the genotype
     * @param referenceBase the reference base
     */
    public VCFGenotypeRecord(String sampleName, Map<String, String> keyValues, List<String> Alleles, PHASE phasing, char referenceBase) {
        mSampleName = sampleName;
        mReferenceBase = referenceBase;
        mFields.putAll(keyValues);
        mAlleleBases.addAll(Alleles);
        phaseType = phasing;
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
        return phaseType;
    }

    public char getReference() {
        return mReferenceBase;
    }

    public List<String> getAllele() {
        return mAlleleBases;
    }

    public Map<String, String> getFields() {
        return mFields;
    }
}
