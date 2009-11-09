package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeCall;

import java.util.*;

/**
 * Split up two call sets into their various concordance sets
 */
public class SNPGenotypeConcordance implements ConcordanceType {

    private double LOD = 5.0;

    private String sample1, sample2;

    public SNPGenotypeConcordance() {}

    public void initialize(Map<String, String> args, Set<String> samples) {
        if ( samples.size() != 2 )
            throw new StingException("SNPGenotype concordance test cannot handle anything other than 2 VCF records");

        if ( args.get("lod") != null )
            LOD = Double.valueOf(args.get("lod"));

        Iterator<String> iter = samples.iterator();
        sample1 = iter.next();
        sample2 = iter.next();
    }

    public String computeConcordance(Map<String, VCFGenotypeCall> samplesToRecords, ReferenceContext ref) {

        VCFGenotypeCall vcfCall1 = samplesToRecords.get(sample1);
        VCFGenotypeCall vcfCall2 = samplesToRecords.get(sample2);
        Variation call1 = (vcfCall1 == null ? null : vcfCall1.toVariation());
        Variation call2 = (vcfCall2 == null ? null : vcfCall2.toVariation());

        // the only reason they would be null is a lack of coverage
        if ( call1 == null || call2 == null ) {
            if ( call1 != null && call1.isSNP() && call1.getNegLog10PError() >= LOD )
                return "set1ConfidentVariantSet2NoCoverage";
            else if ( call2 != null && call2.isSNP() && call2.getNegLog10PError() >= LOD )
                return "set1NoCoverageSet2ConfidentVariant";
            return null;
        }
        if (!(call1 instanceof VariantBackedByGenotype) || !(call2 instanceof VariantBackedByGenotype))
                    throw new StingException("Both ROD tracks must be backed by genotype data. Ensure that your rod(s) contain genotyping information");

        double bestVsRef1 = call1.getNegLog10PError();
        double bestVsRef2 = call2.getNegLog10PError();
        String genotype1 = ((VariantBackedByGenotype)call1).getCalledGenotype().getBases();
        String genotype2 = ((VariantBackedByGenotype)call2).getCalledGenotype().getBases();

        // are they both variant SNPs?
        if ( call1.isSNP() && call2.isSNP() ) {

            // are they both confident calls?
            if ( bestVsRef1 >= LOD && bestVsRef2 >= LOD ) {
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
            else if ( bestVsRef1 < LOD && bestVsRef2 < LOD && bestVsRef1 + bestVsRef2 >= LOD ) {
                return "confidentVariantWhenCombined";
            }

            // only one is confident variant
            else if ( (bestVsRef1 < LOD && bestVsRef2 >= LOD) || (bestVsRef1 >= LOD && bestVsRef2 < LOD) ) {
                return "bothVariantOnlyOneIsConfident";
            }
        }

        // one is variant and the other is ref
        else if ( call1.isSNP() && call2.isReference() && bestVsRef1 >= LOD )
             return "set1VariantSet2Ref";
        else if ( call2.isSNP() && call1.isReference() && bestVsRef2 >= LOD )
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

    public String getInfoName() { return "SG"; }    
}