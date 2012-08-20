package org.broadinstitute.sting.utils.codecs.vcf;

import java.util.Arrays;

/**
 * @author ebanks
 * A class representing a key=value entry for FILTER fields in the VCF header
 */
public class VCFFilterHeaderLine extends VCFSimpleHeaderLine  {

    /**
     * create a VCF filter header line
     *
     * @param name         the name for this header line
     * @param description  the description for this header line
     */
    public VCFFilterHeaderLine(String name, String description) {
        super("FILTER", name, description);
    }

    /**
     * Convenience constructor for FILTER whose description is the name
     * @param name
     */
    public VCFFilterHeaderLine(String name) {
        super("FILTER", name, name);
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     */
    public VCFFilterHeaderLine(String line, VCFHeaderVersion version) {
        super(line, version, "FILTER", Arrays.asList("ID", "Description"));
    }
}