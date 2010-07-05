package org.broad.tribble.vcf;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * @author ebanks
 * A class representing a key=value entry for FILTER fields in the VCF header
 */
public class VCFFilterHeaderLine extends VCFHeaderLine {

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
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     */
    protected VCFFilterHeaderLine(String line, VCFHeaderVersion version) {
        super("FILTER", "", version);
        Map<String,String> mapping = VCFHeaderLineTranslator.parseLine(version,line, Arrays.asList("ID","Description"));
        name = mapping.get("ID");
        description = mapping.get("Description");
    }

    protected String makeStringRep() {
        if (mVersion == VCFHeaderVersion.VCF3_3 || mVersion == VCFHeaderVersion.VCF3_2)
            return String.format("FILTER=%s,\"%s\"", name, description);
        else if (mVersion == VCFHeaderVersion.VCF4_0) {
            Map<String,Object> map = new LinkedHashMap<String,Object>();
            map.put("ID", name);
            map.put("Description", description);
            return "FILTER=" + VCFHeaderLineTranslator.toValue(this.mVersion,map);
        }
        else throw new RuntimeException("Unsupported VCFVersion " + mVersion);
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