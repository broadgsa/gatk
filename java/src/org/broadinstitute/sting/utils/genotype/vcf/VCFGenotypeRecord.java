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
    enum GT_GENOTYPE {
        UNPHASED, PHASED, PHASED_SWITCH_PROB
    }

    // our phasing
    private GT_GENOTYPE phaseType;

    // our reference bases(s)
    private final char mReferenceBase;

    // our bases(s)
    private final List<String> mAlleleBases = new ArrayList<String>();

    // our mapping of the format mFields to values
    private final Map<String, String> mFields = new HashMap<String, String>();
    
    /**
     * create a VCF record
     *
     * @param keyValues     the key values
     * @param Alleles       the alleles, one if we're halpoid, two if we're diploid
     * @param phasing       the phasing of the the genotype
     * @param referenceBase the reference base
     */
    public VCFGenotypeRecord(Map<String, String> keyValues, List<String> Alleles, GT_GENOTYPE phasing, char referenceBase) {
        // validate
        this.mReferenceBase = referenceBase;
        this.mFields.putAll(keyValues);
        this.mAlleleBases.addAll(Alleles);
        this.phaseType = phasing;
    }

    /**
     * determine the phase of the genotype
     *
     * @param phase the string that contains the phase character
     */
    static GT_GENOTYPE determinePhase(String phase) {
        // find the phasing information
        if (phase.equals("/"))
            return GT_GENOTYPE.UNPHASED;
        else if (phase.equals("|"))
            return GT_GENOTYPE.PHASED;
        else if (phase.equals("\\"))
            return GT_GENOTYPE.PHASED_SWITCH_PROB;
        else
            throw new IllegalArgumentException("Unknown genotype phasing parameter");
    }

    /** getter methods */

    public GT_GENOTYPE getPhaseType() {
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
