package org.broadinstitute.sting.utils.genotype.vcf;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;

/**
 * this class writes VCF files
 */
public interface VCFWriter {

    public void writeHeader(VCFHeader header);

    /**
     * attempt to close the VCF file
     */
    public void close();

    public void add(VariantContext vc, byte refBase);
}
