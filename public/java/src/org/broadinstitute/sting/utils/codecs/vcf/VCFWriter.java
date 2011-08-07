package org.broadinstitute.sting.utils.codecs.vcf;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * this class writes VCF files
 */
public interface VCFWriter {

    public void writeHeader(VCFHeader header);

    /**
     * attempt to close the VCF file
     */
    public void close();

    public void add(VariantContext vc);
}