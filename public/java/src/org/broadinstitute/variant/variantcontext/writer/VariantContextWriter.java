package org.broadinstitute.variant.variantcontext.writer;

import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.variantcontext.VariantContext;

/**
 * this class writes VCF files
 */
public interface VariantContextWriter {

    public void writeHeader(VCFHeader header);

    /**
     * attempt to close the VCF file
     */
    public void close();

    public void add(VariantContext vc);
}