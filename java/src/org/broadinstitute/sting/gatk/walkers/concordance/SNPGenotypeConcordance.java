package org.broadinstitute.sting.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Genotype;

import java.util.*;

/**
 * Split up two call sets into their various concordance sets
 */
public class SNPGenotypeConcordance implements ConcordanceType {

    private double Qscore = 30.0;

    private String sample1, sample2;

    public SNPGenotypeConcordance() {}

    public void initialize(Map<String, String> args, Set<String> samples) {
        if ( samples.size() != 2 )
            throw new StingException("SNPGenotype concordance test cannot handle anything other than 2 VCF records");

        if ( args.get("qscore") != null )
            Qscore = Double.valueOf(args.get("qscore"));

        Iterator<String> iter = samples.iterator();
        sample1 = iter.next();
        sample2 = iter.next();
    }

    public String computeConcordance(Map<String, Genotype> samplesToRecords, ReferenceContext ref) {

        Genotype call1 = samplesToRecords.get(sample1);
        Genotype call2 = samplesToRecords.get(sample2);

        // the only reason they would be null is a lack of coverage
        if ( call1 == null || call2 == null ) {
            if ( call1 != null && call1.isPointGenotype() && (10.0 * call1.getNegLog10PError()) >= Qscore )
                return "set1ConfidentVariantSet2NoCoverage";
            else if ( call2 != null && call2.isPointGenotype() && (10.0 * call2.getNegLog10PError()) >= Qscore )
                return "set1NoCoverageSet2ConfidentVariant";
            return null;
        }

        double confidence1 = 10.0 * call1.getNegLog10PError();
        double confidence2 = 10.0 * call2.getNegLog10PError();
        String genotype1 = call1.getBases();
        String genotype2 = call2.getBases();

        // are they both variant SNPs?
        if ( call1.isPointGenotype() && call2.isPointGenotype() ) {

            // are they both confident calls?
            if ( confidence1 >= Qscore && confidence2 >= Qscore ) {
                // same genotype
                if ( genotype1.equals(genotype2) )
                    return "sameConfidentVariant";

                // same allele, different genotype
                else if ( sameVariantAllele(genotype1, genotype2, ref.getBase()) )
                    return "sameVariantAlleleDifferentGenotype";

                // different variant allele
                else
                    return "differentVariantAllele";
            }

            // confident only when combined
            else if ( confidence1 < Qscore && confidence2 < Qscore && confidence1 + confidence2 >= Qscore ) {
                return "confidentVariantWhenCombined";
            }

            // only one is confident variant
            else if ( (confidence1 < Qscore && confidence2 >= Qscore) || (confidence1 >= Qscore && confidence2 < Qscore) ) {
                return "bothVariantOnlyOneIsConfident";
            }
        }

        // one is variant and the other is ref
        else if ( call1.isPointGenotype() && call2.isVariant(ref.getBase()) && confidence1 >= Qscore )
             return "set1VariantSet2Ref";
        else if ( call2.isPointGenotype() && call1.isVariant(ref.getBase()) && confidence2 >= Qscore )
            return "set1RefSet2Variant";

        return null;
    }

    private boolean sameVariantAllele(String genotype1, String genotype2, char ref) {
        if ( genotype1.length() < 2 || genotype2.length() < 2 )
            return genotype1.equals(genotype2);
        char altAllele1 = genotype1.charAt(0) != ref ? genotype1.charAt(0) : genotype1.charAt(1);
        char altAllele2 = genotype2.charAt(0) != ref ? genotype2.charAt(0) : genotype2.charAt(1);
        return altAllele1 == altAllele2;
    }

    public String getInfoName() { return "SnpConcordance"; }    
    public String getInfoDescription() { return getInfoName() + ",1,String,\"SNP concordance test\""; }    
}