package org.broadinstitute.sting.utils.genotype.vcf;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFInfoHeaderLine
 *         <p/>
 *         A class representing a key=value entry for INFO fields in the VCF header
 */
public class VCFInfoHeaderLine extends VCFHeaderLine {

    // the info field types
    public enum INFO_TYPE {
        Integer, Float, String
    }

    private String mName;
    private int mCount;
    private String mDescription;
    private INFO_TYPE mType;


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

    protected String makeStringRep() {
        return String.format("INFO=%s,%d,%s,\"%s\"", mName, mCount, mType.toString(), mDescription);
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
}