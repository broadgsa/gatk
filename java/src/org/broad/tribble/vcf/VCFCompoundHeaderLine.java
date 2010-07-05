package org.broad.tribble.vcf;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * a base class for compound header lines, which include info lines and format lines (so far)
 */
public abstract class VCFCompoundHeaderLine extends VCFHeaderLine {
    public enum SupportedHeaderLineType {
        INFO(true), FORMAT(false);

        public final boolean allowFlagValues;
        SupportedHeaderLineType(boolean flagValues) {
            allowFlagValues = flagValues;
        }

    }

    // the field types
    private String name;
    private int count;
    private String description;
    private VCFHeaderLineType type;

    // access methods
    public String getName() { return name; }
    public int getCount() { return count; }
    public String getDescription() { return description; }
    public VCFHeaderLineType getType() { return type; }

    // our type of line, i.e. format, info, etc
    private final SupportedHeaderLineType lineType;

    // line numerical values are allowed to be unbounded (or unknown), which is
    // marked with a dot (.)
    public static int UNBOUNDED = -1;                       // the value we store internally for unbounded types
    public static String UNBOUNDED_ENCODING_VCF4 = ".";     // the encoding for vcf 4
    public static String UNBOUNDED_ENCODING_VCF3 = "-1";    // the encoding for vcf 3

    /**
     * create a VCF format header line
     *
     * @param name         the name for this header line
     * @param count        the count for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     */
    protected VCFCompoundHeaderLine(String name, int count, VCFHeaderLineType type, String description, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "");
        this.name = name;
        this.count = count;
        this.type = type;
        this.description = description;
        this.lineType = lineType;
    }

    /**
     * create a VCF format header line
     *
     * @param line   the header line
     * @param version      the VCF header version
     *
     */
    protected VCFCompoundHeaderLine(String line, VCFHeaderVersion version, SupportedHeaderLineType lineType) {
        super(lineType.toString(), "", version);
        Map<String,String> mapping = VCFHeaderLineTranslator.parseLine(version,line, Arrays.asList("ID","Number","Type","Description"));
        name = mapping.get("ID");
        count = version == VCFHeaderVersion.VCF4_0 ?
                        mapping.get("Number").equals(UNBOUNDED_ENCODING_VCF4) ? UNBOUNDED : Integer.valueOf(mapping.get("Number")) :
                        mapping.get("Number").equals(UNBOUNDED_ENCODING_VCF3) ? UNBOUNDED : Integer.valueOf(mapping.get("Number"));
        type = VCFHeaderLineType.valueOf(mapping.get("Type"));
        if (type == VCFHeaderLineType.Flag && !allowFlagValues())
            throw new IllegalArgumentException("Flag is an unsupported type for this kind of field");
        description = mapping.get("Description");
        this.lineType = lineType;
    }

    /**
     * make a string representation of this header line
     * @return a string representation
     */
    protected String makeStringRep() {
        if (mVersion == VCFHeaderVersion.VCF3_3 || mVersion == VCFHeaderVersion.VCF3_2)
            return String.format(lineType.toString()+"=%s,%d,%s,\"%s\"", name, count, type.toString(), description);
        else if (mVersion == VCFHeaderVersion.VCF4_0) {
            Map<String,Object> map = new LinkedHashMap<String,Object>();
            map.put("ID", name);
            // TODO: this next line should change when we have more than two used encoding schemes
            map.put("Number", count == UNBOUNDED ? (mVersion == VCFHeaderVersion.VCF4_0 ? UNBOUNDED_ENCODING_VCF4 : UNBOUNDED_ENCODING_VCF3) : count);
            map.put("Type", type);
            map.put("Description", description);
            return lineType.toString() + "=" + VCFHeaderLineTranslator.toValue(this.mVersion,map);
        }
        else throw new RuntimeException("Unsupported VCFVersion " + mVersion);
    }

    /**
     * returns true if we're equal to another compounder header line
     * @param o a compound header line
     * @return true if equal
     */
    public boolean equals(Object o) {
        if ( !(o instanceof VCFCompoundHeaderLine) )
            return false;
        VCFCompoundHeaderLine other = (VCFCompoundHeaderLine)o;
        return name.equals(other.name) &&
               count == other.count &&
               description.equals(other.description) &&
               type == other.type;
    }

    /**
     * do we allow flag (boolean) values? (i.e. booleans where you don't have specify the value, AQ means AQ=true)
     * @return true if we do, false otherwise
     */
    abstract boolean allowFlagValues();

}
