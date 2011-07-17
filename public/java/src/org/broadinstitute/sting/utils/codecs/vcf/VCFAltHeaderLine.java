package org.broadinstitute.sting.utils.codecs.vcf;

/**
 * @author ebanks
 * A class representing a key=value entry for ALT fields in the VCF header
 */
public class VCFAltHeaderLine extends VCFSimpleHeaderLine  {

    /**
     * create a VCF filter header line
     *
     * @param name         the name for this header line
     * @param description  the description for this header line
     */
    public VCFAltHeaderLine(String name, String description) {
        super(name, description, SupportedHeaderLineType.ALT);
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     */
    protected VCFAltHeaderLine(String line, VCFHeaderVersion version) {
        super(line, version, SupportedHeaderLineType.ALT);
    }
}