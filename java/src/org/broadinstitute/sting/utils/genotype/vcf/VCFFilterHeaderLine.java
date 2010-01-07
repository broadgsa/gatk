package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.Utils;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFFilterHeaderLine
 *         <p/>
 *         A class representing a key=value entry for FILTER fields in the VCF header
 */
public class VCFFilterHeaderLine extends VCFHeaderLine {

    private String mName;
    private String mDescription;


    /**
     * create a VCF filter header line
     *
     * @param name         the name for this header line
     * @param description  the description for this header line
     */
    public VCFFilterHeaderLine(String name, String description) {
        super("FILTER", "");
        mName = name;
        mDescription = description;
    }

    /**
     * create a VCF info header line
     *
     * @param line   the header line
     */
    protected VCFFilterHeaderLine(String line) {
        super("FILTER", "");
        String[] pieces = line.split(",");
        if ( pieces.length < 2 )
            throw new IllegalArgumentException("There are too few values in the VCF FILTER header line: " + line);

        mName = pieces[0];
        mDescription = Utils.trim(pieces[1], '"');
        // just in case there were some commas in the description
        for (int i = 1; i < pieces.length; i++)
            mDescription += "," + Utils.trim(pieces[i], '"');
    }

    protected String makeStringRep() {
        return String.format("FILTER=%s,\"%s\"", mName, mDescription);
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFFilterHeaderLine) )
            return false;
        VCFFilterHeaderLine other = (VCFFilterHeaderLine)o;
        return mName.equals(other.mName) &&
               mDescription.equals(other.mDescription);
    }
}