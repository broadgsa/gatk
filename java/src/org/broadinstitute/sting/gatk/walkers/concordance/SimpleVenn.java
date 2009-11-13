package org.broadinstitute.sting.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeCall;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;

/**
 * Split up two call sets into their various "Venn diagram" sets
 */
public class SimpleVenn implements ConcordanceType {

    private String sample1, sample2;

    public SimpleVenn() {}

    public void initialize(Map<String, String> args, Set<String> samples) {
        if ( samples.size() != 2 )
            throw new StingException("SimpleVenn concordance test cannot handle anything other than 2 VCF records");        

        Iterator<String> iter = samples.iterator();
        sample1 = iter.next();
        sample2 = iter.next();
    }

    public String computeConcordance(Map<String, VCFGenotypeCall> samplesToRecords, ReferenceContext ref) {

        VCFGenotypeCall call1 = samplesToRecords.get(sample1);
        VCFGenotypeCall call2 = samplesToRecords.get(sample2);

        if ( call1 == null && call2 == null )
            return null;

        // set 1 only
        if ( call2 == null )
            return sample1 + "_only";

        // set 2 only
        else if ( call1 == null )
            return sample2 + "_only";

        // at this point we know that neither is null, so now we need to test for alternate allele concordance
        Variation callV1 = call1.toVariation(ref.getBase());
        Variation callV2 = call2.toVariation(ref.getBase());

        // we can't really deal with multi-allelic variants
        if ( callV1.isBiallelic() && callV2.isBiallelic() ) {
            // intersection (concordant)
            if ( callV1.getAlternativeBaseForSNP() == callV2.getAlternativeBaseForSNP() )
                return "concordant";
            // intersection (discordant)
            else
                return "discordant";
        }

        return "concordant";
    }

    public String getInfoName() { return "Venn"; }
}
