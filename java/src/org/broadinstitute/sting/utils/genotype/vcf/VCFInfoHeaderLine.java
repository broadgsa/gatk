package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.Utils;


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
        Integer, Float, String, Character, Flag
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

    /**
     * create a VCF info header line
     *
     * @param line   the header line
     */
    protected VCFInfoHeaderLine(String line) {
        super("INFO", "");
        String[] pieces = line.split(",");
        if ( pieces.length < 4 )
            throw new IllegalArgumentException("There are too few values in the VCF INFO header line: " + line);

        mName = pieces[0];
        mCount = Integer.valueOf(pieces[1]);
        mType = INFO_TYPE.valueOf(pieces[2]);
        mDescription = Utils.trim(pieces[3], '"');
        // just in case there were some commas in the description
        for (int i = 4; i < pieces.length; i++)
            mDescription += "," + Utils.trim(pieces[i], '"');
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