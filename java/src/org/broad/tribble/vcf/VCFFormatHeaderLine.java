package org.broad.tribble.vcf;

import org.broad.tribble.util.ParsingUtils;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFFormatHeaderLine
 *         <p/>
 *         A class representing a key=value entry for genotype FORMAT fields in the VCF header
 */
public class VCFFormatHeaderLine extends VCFHeaderLine {

    // the format field types
    public enum FORMAT_TYPE {
        Integer, Float, String;
        public Object convert(String value) {
            switch (this) {
                case Integer:
                    // Take care of case when string might have say "43.0" but we request an Integer.
                    return Math.round(java.lang.Float.valueOf(value)); // the java.lang is needed since we use Integer as a enum name
                case Float:
                    return java.lang.Float.valueOf(value);
                case String:
                    return value;
                default:
                    throw new IllegalStateException("field." + this + " doesn't have a set conversion approach");
            }
        }
    }

    private String mName;
    private int mCount;
    private String mDescription;
    private FORMAT_TYPE mType;

    // info line numerical values are allowed to be unbounded (or unknown), which is
    // marked with a dot (.)
    public static int UNBOUNDED = -1;
    public static String UNBOUNDED_ENCODING_VCF4 = ".";
    public static String UNBOUNDED_ENCODING_VCF3 = "-1";

    /**
     * create a VCF format header line
     *
     * @param name         the name for this header line
     * @param count        the count for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     */
    public VCFFormatHeaderLine(String name, int count, FORMAT_TYPE type, String description) {
        super("FORMAT", "");
        mName = name;
        mCount = count;
        mType = type;
        mDescription = description;
    }

    /**
     * create a VCF format header line
     *
     * @param line   the header line
     * @param version      the VCF header version
     *
     */
    protected VCFFormatHeaderLine(String line, VCFHeaderVersion version) {
        super("FORMAT", "", version);
        Map<String,String> mapping = VCFHeaderLineTranslator.parseLine(version,line, Arrays.asList("ID","Number","Type","Description"));
        mName = mapping.get("ID");
        mCount = version == VCFHeaderVersion.VCF4_0 ?
                        mapping.get("Number").equals(UNBOUNDED_ENCODING_VCF4) ? UNBOUNDED : Integer.valueOf(mapping.get("Number")) :
                        mapping.get("Number").equals(UNBOUNDED_ENCODING_VCF3) ? UNBOUNDED : Integer.valueOf(mapping.get("Number"));
        mType = FORMAT_TYPE.valueOf(mapping.get("Type"));
        mDescription = mapping.get("Description");
    }

    protected String makeStringRep() {
        if (mVersion == VCFHeaderVersion.VCF3_3 || mVersion == VCFHeaderVersion.VCF3_2)
            return String.format("FORMAT=%s,%d,%s,\"%s\"", mName, mCount, mType.toString(), mDescription);
        else if (mVersion == VCFHeaderVersion.VCF4_0) {
            Map<String,Object> map = new LinkedHashMap<String,Object>();
            map.put("ID",mName);
            map.put("Number",mCount == UNBOUNDED ? (mVersion == VCFHeaderVersion.VCF4_0 ? UNBOUNDED_ENCODING_VCF4 : UNBOUNDED_ENCODING_VCF3) : mCount);
            map.put("Type",mType);
            map.put("Description",mDescription);
            return "FORMAT=" + VCFHeaderLineTranslator.toValue(this.mVersion,map);
        }
        else throw new RuntimeException("Unsupported VCFVersion " + mVersion);
    }

    public String getName() { return mName; }
    public int getCount() { return mCount; }
    public String getDescription() { return mDescription; }
    public FORMAT_TYPE getType() { return mType; }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFFormatHeaderLine) )
            return false;
        VCFFormatHeaderLine other = (VCFFormatHeaderLine)o;
        return mName.equals(other.mName) &&
               mCount == other.mCount &&
               mDescription.equals(other.mDescription) &&
               mType == other.mType;
    }

    public String getmName() {
        return mName;
    }

    public int getmCount() {
        return mCount;
    }

    public String getmDescription() {
        return mDescription;
    }

    public FORMAT_TYPE getmType() {
        return mType;
    }
}