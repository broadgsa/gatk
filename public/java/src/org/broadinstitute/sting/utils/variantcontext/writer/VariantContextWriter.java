package org.broadinstitute.sting.utils.variantcontext.writer;

import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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