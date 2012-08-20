package org.broadinstitute.sting.utils.codecs.vcf;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


/**
 * @author ebanks
 * A class representing a key=value entry for simple VCF header types
 */
public class VCFSimpleHeaderLine extends VCFHeaderLine implements VCFIDHeaderLine {

    private String name;
    private Map<String, String> genericFields = new LinkedHashMap<String, String>();

    /**
     * create a VCF filter header line
     *
     * @param key            the key for this header line
     * @param name           the name for this header line
     * @param description    description for this header line
     */
    public VCFSimpleHeaderLine(String key, String name, String description) {
        super(key, "");
        Map<String, String> map = new LinkedHashMap<String, String>(1);
        map.put("Description", description);
        initialize(name, map);
    }

    /**
     * create a VCF info header line
     *
     * @param line      the header line
     * @param version   the vcf header version
     * @param key            the key for this header line
     * @param expectedTagOrdering the tag ordering expected for this header line
     */
    public VCFSimpleHeaderLine(final String line, final VCFHeaderVersion version, final String key, final List<String> expectedTagOrdering) {
        this(key, VCFHeaderLineTranslator.parseLine(version, line, expectedTagOrdering), expectedTagOrdering);
    }

    public VCFSimpleHeaderLine(final String key, final Map<String, String> mapping, final List<String> expectedTagOrdering) {
        super(key, "");
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
        return getKey() + "=" + VCFHeaderLine.toStringEncoding(map);
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
}