package org.broadinstitute.sting.utils.genotype.vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


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
    private final char reference;

    // our bases(s)
    private final List<String> bases = new ArrayList<String>();

    // our mapping of the format fields to values
    private final Map<String, String> fields = new HashMap<String, String>();

    // our pattern matching for the genotype fields
    private static final Pattern basicSplit = Pattern.compile("([0-9]*)([\\\\|\\/])([0-9]*):(\\S*)");
    
    /**
     * generate a VCF genotype record, given it's format string, the genotype string, and allele info
     *
     * @param formatString   the format string for this record, which contains the keys for the genotype parameters
     * @param genotypeString contains the phasing information, allele information, and values for genotype parameters
     * @param altAlleles     the alternate allele string array, which we index into based on the field parameters
     * @param referenceBase  the reference base
     */
    protected VCFGenotypeRecord(String formatString, String genotypeString, String altAlleles[], char referenceBase) {
        reference = referenceBase;
        // check that the first format field is GT, which is required
        String keys[] = formatString.split(":");
        if (keys.length < 0 || !keys[0].equals("GT"))
            throw new IllegalArgumentException("The format string must have fields, and the first must be GT (genotype)");

        // find the values for each of the keys, of which the GT field should be the first
        Matcher match = basicSplit.matcher(genotypeString);
        if (!match.matches() || match.groupCount() < 3)
            throw new IllegalArgumentException("Unable to match genotype string to expected regex");

        // add the alternate base (which can be ref by specifying 0)
        addAllele(match.group(1), altAlleles, referenceBase);

        determinePhase(match.group(2));

        // do we have a second alt allele?
        if (match.group(3).length() > 0) {
            addAllele(match.group(3), altAlleles, referenceBase);
        }

        // check to see what other records we have
        if (match.groupCount() == 4) {
            // make sure we'll have enough occurances
            String tokens[] = match.group(4).split(":{1}"); // the {1} was required, since string.split does a greedy match of the specified regex, like :+
            int keyIndex = 1;
            for (String token: tokens) {
                this.fields.put(keys[keyIndex],token);
                keyIndex++;
            }
            if (keyIndex + 1 == tokens.length) fields.put(keys[++keyIndex],""); // if the last value is blank, split will leave it off
            if (keyIndex == 1 && match.group(4).contains(":")) {
                // there was a string of all semicolons, split doesn't handle this well (or at all)
                while(keyIndex < keys.length) this.fields.put(keys[keyIndex++],"");
            }
        }

    }

    /**
     * add an alternate allele to the list of alleles we have
     *
     * @param alleleNumber  the allele number, as a string
     * @param altAlleles    the list of alternate alleles
     * @param referenceBase the reference base
     */
    private void addAllele(String alleleNumber, String[] altAlleles, char referenceBase) {
        if (Integer.valueOf(alleleNumber) == 0)
            bases.add(String.valueOf(referenceBase));
        else
            bases.add(altAlleles[Integer.valueOf(alleleNumber) - 1]);
    }

    /**
     * determine the phase of the genotype
     *
     * @param phase the string that contains the phase character
     */
    private void determinePhase(String phase) {
        // find the phasing information
        if (phase.equals("/"))
            phaseType = GT_GENOTYPE.UNPHASED;
        else if (phase.equals("|"))
            phaseType = GT_GENOTYPE.PHASED;
        else if (phase.equals("\\"))
            phaseType = GT_GENOTYPE.PHASED_SWITCH_PROB;
        else
            throw new IllegalArgumentException("Unknown genotype phasing parameter");
    }

    /** getter methods */

    public GT_GENOTYPE getPhaseType() {
        return phaseType;
    }

    public char getReference() {
        return reference;
    }

    public List<String> getAllele() {
        return bases;
    }

    public Map<String, String> getFields() {
        return fields;
    }
}
