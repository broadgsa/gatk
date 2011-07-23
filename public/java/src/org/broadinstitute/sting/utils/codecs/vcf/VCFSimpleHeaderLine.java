package org.broadinstitute.sting.utils.codecs.vcf;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * @author ebanks
 * A class representing a key=value entry for simple VCF header types
 */
public abstract class VCFSimpleHeaderLine extends VCFHeaderLine implements VCFNamedHeaderLine  {

    public enum SupportedHeaderLineType {
        FILTER, ALT;
    }

    private String name;
    private String description;

    // our type of line, i.e. filter, alt, etc
    private final SupportedHeaderLineType lineType;


    /**
     * create a VCF filter header line
     *
     * @param name         the name for this header line
     * @param description  the description for this header line
     * @param lineType     the header line type
     */
    public VCFSimpleHeaderLine(String name, String description, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.lineType = lineType;
        this.name = name;
        this.description = description;

        if ( name == null || description == null )
            throw new IllegalArgumentException(String.format("Invalid VCFSimpleHeaderLine: key=%s name=%s desc=%s", super.getKey(), name, description ));
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     * @param lineType     the header line type
     */
    protected VCFSimpleHeaderLine(String line, VCFHeaderVersion version, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.lineType = lineType;
        Map<String,String> mapping = VCFHeaderLineTranslator.parseLine(version,line, Arrays.asList("ID","Description"));
        name = mapping.get("ID");
        description = mapping.get("Description");
        if ( description == null && ALLOW_UNBOUND_DESCRIPTIONS ) // handle the case where there's no description provided 
            description = UNBOUND_DESCRIPTION;
    }

    protected String toStringEncoding() {
        Map<String,Object> map = new LinkedHashMap<String,Object>();
        map.put("ID", name);
        map.put("Description", description);
        return lineType.toString() + "=" + VCFHeaderLine.toStringEncoding(map);
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFSimpleHeaderLine) )
            return false;
        VCFSimpleHeaderLine other = (VCFSimpleHeaderLine)o;
        return name.equals(other.name) &&
               description.equals(other.description);
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }
}