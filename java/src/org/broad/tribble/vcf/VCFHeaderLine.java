package org.broad.tribble.vcf;



/**
 * @author ebanks
 *         <p/>
 *         Class VCFHeaderLine
 *         <p/>
 *         A class representing a key=value entry in the VCF header
 */
public class VCFHeaderLine implements Comparable {

    private String stringRep = null;
    private String mKey = null;
    private String mValue = null;
    protected VCFHeaderVersion mVersion = null;

    /**
     * create a VCF header line
     *
     * @param key     the key for this header line
     * @param value   the value for this header line
     */
    public VCFHeaderLine(String key, String value, VCFHeaderVersion version) {
        mKey = key;
        mValue = value;
        mVersion = version;
    }

    /**
     * create a VCF header line
     *
     * @param key     the key for this header line
     * @param value   the value for this header line
     */
    public VCFHeaderLine(String key, String value) {
        mKey = key;
        mValue = value;
        mVersion = VCFHeaderVersion.VCF3_3;
    }

    /**
     * Get the key
     *
     * @return the key
     */
    public String getKey() {
        return mKey;
    }

    /**
     * Set the key
     *
     * @param key     the key for this header line
     */
    public void setKey(String key) {
        mKey = key;
        stringRep = null;
    }

    /**
     * Get the value
     *
     * @return the value
     */
    public String getValue() {
        return mValue;
    }

    /**
     * Set the value
     *
     * @param value     the value for this header line
     */
    public void setValue(String value) {
        mValue = value;
        stringRep = null;
    }

    public String toString() {
        if ( stringRep == null )
            stringRep = makeStringRep();
        return stringRep;
    }

    protected String makeStringRep() {
        return mKey + "=" + mValue;
    }

    public boolean equals(Object o) {
        if ( !(o instanceof VCFHeaderLine) )
            return false;
        return mKey.equals(((VCFHeaderLine)o).getKey()) && mValue.equals(((VCFHeaderLine)o).getValue());
    }

    public int compareTo(Object other) {
        return toString().compareTo(other.toString());
    }
}