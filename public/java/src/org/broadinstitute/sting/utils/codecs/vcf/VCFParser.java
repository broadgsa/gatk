package org.broadinstitute.sting.utils.codecs.vcf;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;

import java.util.List;


/**
 * All VCF codecs need to implement this interface so that we can perform lazy loading.
 */
public interface VCFParser {

    /**
     * create a genotype map
     * @param str the string
     * @param alleles the list of alleles
     * @param chr chrom
     * @param pos position
     * @return a mapping of sample name to genotype object
     */
    public GenotypesContext createGenotypeMap(String str, List<Allele> alleles, String chr, int pos);

}
