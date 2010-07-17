package org.broad.tribble.vcf;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * @author ebanks
 * A class representing a key=value entry for FILTER fields in the VCF header
 */
public class VCFFilterHeaderLine extends VCFHeaderLine implements VCFNamedHeaderLine  {

    private String name;
    private String description;


    /**
     * create a VCF filter header line
     *
     * @param name         the name for this header line
     * @param description  the description for this header line
     */
    public VCFFilterHeaderLine(String name, String description) {
        super("FILTER", "");
        this.name = name;
        this.description = description;

        if ( name == null || description == null )
            throw new IllegalArgumentException(String.format("Invalid VCFCompoundHeaderLine: key=%s name=%s desc=%s", super.getKey(), name, description ));
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     */
    protected VCFFilterHeaderLine(String line, VCFHeaderVersion version) {
        super("FILTER", "");
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
        return "FILTER=" + VCFHeaderLine.toStringEncoding(map);
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFFilterHeaderLine) )
            return false;
        VCFFilterHeaderLine other = (VCFFilterHeaderLine)o;
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