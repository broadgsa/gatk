package org.broadinstitute.sting.gatk.refdata.features.vcf4;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

/**
 * simple variant context wrapped as VCF4
 */
public class VCF4Record extends VariantContext implements Feature {
    /**
     * create a VCF4Record, which is really a variant context
     * @param name the name of the record
     * @param loc it's location
     * @param alleles the set of alleles
     * @param genotypes any genotypes for this record
     * @param negLog10PError the probability of being a wrong call
     * @param filters the set of filters applied to this variant
     * @param attributes any other attributes
     */
    public VCF4Record(String name, GenomeLoc loc, Collection<Allele> alleles, Map<String, Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        super(name, loc, alleles, genotypes, negLog10PError, filters, attributes);
    }

    @Override
    public String getChr() {
        return getLocation().getContig();
    }

    @Override
    public int getStart() {
        return (int)getLocation().getStart();
    }

    @Override
    public int getEnd() {
        return (int)getLocation().getStop();
    }
}
