package org.broadinstitute.sting.utils.codecs.vcf;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFFormatHeaderLine
 *         <p/>
 *         A class representing a key=value entry for genotype FORMAT fields in the VCF header
 */
public class VCFFormatHeaderLine extends VCFCompoundHeaderLine {

    public VCFFormatHeaderLine(String name, int count, VCFHeaderLineType type, String description) {
        super(name, count, type, description, SupportedHeaderLineType.FORMAT);
        if (type == VCFHeaderLineType.Flag)
            throw new IllegalArgumentException("Flag is an unsupported type for format fields");
    }

    public VCFFormatHeaderLine(String name, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
        super(name, count, type, description, SupportedHeaderLineType.FORMAT);
    }

    protected VCFFormatHeaderLine(String line, VCFHeaderVersion version) {
        super(line, version, SupportedHeaderLineType.FORMAT);
    }

    // format fields do not allow flag values (that wouldn't make much sense, how would you encode this in the genotype).
    @Override
    boolean allowFlagValues() {
        return false;
    }
}