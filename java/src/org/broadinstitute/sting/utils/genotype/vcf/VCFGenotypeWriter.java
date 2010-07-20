package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;

import java.io.File;

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
     * @param header  the header
     */
    public void writeHeader(VCFHeader header);

    /**
     * Add a given VCF file to the writer.
     * @param file  file from which to add records
     */
    public void append(File file);

}
