package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.Utils;


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
        Integer, Float, String
    }

    private String mName;
    private int mCount;
    private String mDescription;
    private FORMAT_TYPE mType;


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
     */
    protected VCFFormatHeaderLine(String line) {
        super("FORMAT", "");
        String[] pieces = line.split(",");
        if ( pieces.length < 4 )
            throw new IllegalArgumentException("There are too few values in the VCF FORMAT header line: " + line);

        mName = pieces[0];
        mCount = Integer.valueOf(pieces[1]);
        mType = FORMAT_TYPE.valueOf(pieces[2]);
        mDescription = Utils.trim(pieces[3], '"');
        // just in case there were some commas in the description
        for (int i = 4; i < pieces.length; i++)
            mDescription += "," + Utils.trim(pieces[i], '"');
    }

    protected String makeStringRep() {
        return String.format("FORMAT=%s,%d,%s,\"%s\"", mName, mCount, mType.toString(), mDescription);
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
}