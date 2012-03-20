package org.broadinstitute.sting.utils.codecs.vcf;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


/**
 * @author ebanks
 * A class representing a key=value entry for simple VCF header types
 */
public class VCFSimpleHeaderLine extends VCFHeaderLine implements VCFIDHeaderLine {

    public enum SupportedHeaderLineType {
        FILTER, GENERIC;
    }

    private String name;
    private Map<String, String> genericFields = new LinkedHashMap<String, String>();

    
    // our type of line, i.e. filter, alt, etc
    private final SupportedHeaderLineType lineType;


    /**
     * create a VCF filter header line
     *
     * @param name           the name for this header line
     * @param genericFields  other fields for this header line
     * @param lineType       the header line type
     */
    public VCFSimpleHeaderLine(String name, Map<String, String> genericFields, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.lineType = lineType;
        initialize(name, genericFields);
    }

    /**
     * create a VCF filter header line
     *
     * @param name           the name for this header line
     * @param description    description for this header line
     * @param lineType       the header line type
     */
    public VCFSimpleHeaderLine(String name, String description, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.lineType = lineType;
        Map<String, String> map = new LinkedHashMap<String, String>(1);
        map.put("Description", description);
        initialize(name, map);
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     * @param lineType  the header line type
     * @param expectedTagOrdering the tag ordering expected for this header line
     */
    protected VCFSimpleHeaderLine(String line, VCFHeaderVersion version, SupportedHeaderLineType lineType, List<String> expectedTagOrdering) {
        super(lineType.toString(), "");
        this.lineType = lineType;
        Map<String, String> mapping = VCFHeaderLineTranslator.parseLine(version, line, expectedTagOrdering);
        name = mapping.get("ID");
        initialize(name, mapping);
    }

    protected void initialize(String name, Map<String, String> genericFields) {
        if ( name == null || genericFields == null || genericFields.isEmpty() )
            throw new IllegalArgumentException(String.format("Invalid VCFSimpleHeaderLine: key=%s name=%s", super.getKey(), name));

        this.name = name;
        this.genericFields.putAll(genericFields);
    }

    protected String toStringEncoding() {
        Map<String, Object> map = new LinkedHashMap<String, Object>();
        map.put("ID", name);
        map.putAll(genericFields);
        return lineType.toString() + "=" + VCFHeaderLine.toStringEncoding(map);
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFSimpleHeaderLine) )
            return false;
        VCFSimpleHeaderLine other = (VCFSimpleHeaderLine)o;
        if ( !name.equals(other.name) || genericFields.size() != other.genericFields.size() )
            return false;
        for ( Map.Entry<String, String> entry : genericFields.entrySet() ) {
            if ( !entry.getValue().equals(other.genericFields.get(entry.getKey())) )
                return false;
        }
        
        return true;       
    }

    public String getID() {
        return name;
    }

    public Map<String, String> getGenericFields() {
        return genericFields;
    }
}