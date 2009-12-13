package org.broadinstitute.sting.utils.genotype.vcf;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFFormatHeaderLine
 *         <p/>
 *         A class representing a key=value entry for genotype FORMAT fields in the VCF header
 */
public class VCFFormatHeaderLine extends VCFHeaderLine {

    // the info field types
    public enum INFO_TYPE {
        Integer, Float, String
    }

    private String mName;
    private int mCount;
    private String mDescription;
    private INFO_TYPE mType;


    /**
     * create a VCF format header line
     *
     * @param name         the name for this header line
     * @param count        the count for this header line
     * @param type         the type for this header line
     * @param description  the description for this header line
     */
    public VCFFormatHeaderLine(String name, int count, INFO_TYPE type, String description) {
        super("FORMAT", "");
        mName = name;
        mCount = count;
        mType = type;
        mDescription = description;
    }

    protected String makeStringRep() {
        return String.format("FORMAT=%s,%d,%s,\"%s\"", mName, mCount, mType.toString(), mDescription);
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFFormatHeaderLine) )
            return false;
        VCFFormatHeaderLine other = (VCFFormatHeaderLine)o;
        return mName.equals(other.mName) &&
               mCount == other.mCount &&
               mDescription.equals(other.mDescription) &&
               mType == other.mType;
    }
}