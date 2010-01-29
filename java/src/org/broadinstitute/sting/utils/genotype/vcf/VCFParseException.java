package org.broadinstitute.sting.utils.genotype.vcf;

/**
 * an exception to funnel all parsing exceptions into; this way we can emit the line we choked on as well
 */
public class VCFParseException extends RuntimeException {
    public VCFParseException(String message) {
        super(message);
    }

    public VCFParseException(String message, Throwable cause) {
        super(message, cause);
    }
}
