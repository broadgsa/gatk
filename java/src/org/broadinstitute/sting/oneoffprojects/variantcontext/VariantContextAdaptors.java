package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.gatk.refdata.rodDbSNP;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;


public class VariantContextAdaptors {
    public static VariantContext dbsnp2VariantContext(rodDbSNP dbsnp) {
        VariantContext vc = new VariantContext(dbsnp.getLocation());

        // add the reference allele
        Allele refAllele = new Allele(dbsnp.getReference(), true);
        vc.addAllele(refAllele);

        // add all of the alt alleles
        for ( String alt : dbsnp.getAlternateAlleleList() )
            vc.addAllele(new Allele(alt, false));

        return vc;
    }
}