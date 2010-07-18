package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;

import java.util.Set;

/**
 * An extension of the GenotypeWriter interface with support
 * for adding header lines.
 *
 * @author mhanna
 * @version 0.1
 */
public interface VCFGenotypeWriter extends GenotypeWriter {
    /**
     * initialize this VCF header
     *
     * @param sampleNames  the sample names
     * @param headerInfo  the optional header fields
     */
    public void writeHeader(Set<String> sampleNames, Set<VCFHeaderLine> headerInfo);

    /**
     * Add a given VCF record to the given output.
     * @param vcfRecord Record to add.
     */
    public void addRecord(VCFRecord vcfRecord);
}
