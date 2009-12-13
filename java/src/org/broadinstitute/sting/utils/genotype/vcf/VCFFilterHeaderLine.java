package org.broadinstitute.sting.utils.genotype.vcf;


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