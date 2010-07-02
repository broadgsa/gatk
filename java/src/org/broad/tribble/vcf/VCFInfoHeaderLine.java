package org.broad.tribble.vcf;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFInfoHeaderLine
 *         <p/>
 *         A class representing a key=value entry for INFO fields in the VCF header
 */
public class VCFInfoHeaderLine extends VCFHeaderLine implements VCFNamedHeaderLine {

    // the info field types
    public enum INFO_TYPE {
        Integer, Float, String, Character, Flag;

        public Object convert(String value) {
            switch (this) {
                case Integer:
                    return java.lang.Integer.valueOf(value); // the java.lang is needed since we use Integer as a enum name
                case Float:
                    return java.lang.Float.valueOf(value);
                case String:
                case Character:
                    return value;
                case Flag:
                    return value.equals("0") ? false : true;
                default:
                    throw new IllegalStateException("INFO_TYPE." + this + " doesn't have a set conversion approach");
            }
        }
    }

    private String mName;
    private int mCount;
    private String mDescription;
    private INFO_TYPE mType;


    // info line numerical values are allowed to be unbounded (or unknown), which is
    // marked with a dot (.)
    public static int UNBOUNDED = -1;
    public static String UNBOUNDED_ENCODING_VCF4 = ".";
    public static String UNBOUNDED_ENCODING_VCF3 = "-1";

    /**
     * create a VCF info header line
     *
     * @param name         the name for this header line
     * @param count        the count for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     */
    public VCFInfoHeaderLine(String name, int count, INFO_TYPE type, String description) {
        super("INFO", "");
        mName = name;
        mCount = count;
        mType = type;
        mDescription = description;
    }

    /**
     * create a VCF info header line
     *
     * @param line   the header line
     * @param version the VCF version
     */
    protected VCFInfoHeaderLine(String line, VCFHeaderVersion version) {
        super("INFO", "", version);
        Map<String,String> mapping = VCFHeaderLineTranslator.parseLine(version,line, Arrays.asList("ID","Number","Type","Description"));
        mName = mapping.get("ID");
        mCount = version == VCFHeaderVersion.VCF4_0 ?
                mapping.get("Number").equals(UNBOUNDED_ENCODING_VCF4) ? UNBOUNDED : Integer.valueOf(mapping.get("Number")) :
                mapping.get("Number").equals(UNBOUNDED_ENCODING_VCF3) ? UNBOUNDED : Integer.valueOf(mapping.get("Number"));
        mType = INFO_TYPE.valueOf(mapping.get("Type"));
        mDescription = mapping.get("Description");
    }

    protected String makeStringRep() {
        if (mVersion == VCFHeaderVersion.VCF3_3 || mVersion == VCFHeaderVersion.VCF3_2)
            return String.format("INFO=%s,%d,%s,\"%s\"", mName, mCount, mType.toString(), mDescription);
        else if (mVersion == VCFHeaderVersion.VCF4_0) {
            Map<String,Object> map = new LinkedHashMap<String,Object>();
            map.put("ID",mName);
            map.put("Number",mCount == UNBOUNDED ? (mVersion == VCFHeaderVersion.VCF4_0 ? UNBOUNDED_ENCODING_VCF4 : UNBOUNDED_ENCODING_VCF3) : mCount);
            map.put("Type",mType);
            map.put("Description",mDescription);
            return "INFO=" + VCFHeaderLineTranslator.toValue(this.mVersion,map);
        }
        else throw new RuntimeException("Unsupported VCFVersion " + mVersion);
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFInfoHeaderLine) )
            return false;
        VCFInfoHeaderLine other = (VCFInfoHeaderLine)o;
        return mName.equals(other.mName) &&
               mCount == other.mCount &&
               mDescription.equals(other.mDescription) &&
               mType == other.mType;
    }

    @Override
    public String getmName() {
        return mName;
    }

    public int getmCount() {
        return mCount;
    }

    public String getmDescription() {
        return mDescription;
    }

    public INFO_TYPE getmType() {
        return mType;
    }
}