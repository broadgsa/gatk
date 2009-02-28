/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.sam;


import java.util.*;

/**
 * Header information about a reference sequence.
 */
public class SAMSequenceRecord
{
    private String mSequenceName = null;
    private int mSequenceIndex = -1;
    private int mSequenceLength = 0;
    private Map<String, Object> mAttributes = null;
    public static final String SEQUENCE_NAME_TAG = "SN";
    public static final String SEQUENCE_LENGTH_TAG = "LN";
    public static final String MD5_TAG = "M5";
    public static final String ASSEMBLY_TAG = "AS";
    public static final String URI_TAG = "UR";
    public static final String SPECIES_TAG = "SP";

    public SAMSequenceRecord(final String name) {
        mSequenceName = name;
    }

    public String getSequenceName() {
        return mSequenceName;
    }

    public int getSequenceLength() {
        return mSequenceLength;
    }

    public void setSequenceLength(final int value) {
        mSequenceLength = value;
    }

    public String getAssembly() {
        return (String) getAttribute("AS");
    }

    public void setAssembly(final String value) {
        setAttribute("AS", value);
    }

    public String getSpecies() {
        return (String) getAttribute("SP");
    }

    public void setSpecies(final String value) {
        setAttribute("SP", value);
    }

    public Object getAttribute(final String key) {
        if (mAttributes == null) {
            return null;
        }
        return mAttributes.get(key);
    }

    public void setAttribute(final String key, final Object value) {
        if (mAttributes == null) {
            mAttributes = new HashMap<String, Object>();
        }
        mAttributes.put(key, value);
    }

    public Set<Map.Entry<String, Object>> getAttributes() {
        if (mAttributes == null) {
            return null;
        }
        return mAttributes.entrySet();
    }

    // Private state used only by SAM implementation.
    int getSequenceIndex() {
        return mSequenceIndex;
    }

    // Private state used only by SAM implementation.
    void setSequenceIndex(final int value) {
        mSequenceIndex = value;
    }

    /**
     * Looser comparison than equals().  If one SAMSequenceRecord has an attribute that the other does not
     * have, that is not considered inequality.  However, if they both have an attribute, but have different
     * values for that atttribute, then they are considered unequal.  This results in an intransitive equality test,
     * i.e. a.isSameSequence(b) && b.isSameSequence(c) does not necessarily imply a.isSameSequence(c)
     */
    public boolean isSameSequence(final SAMSequenceRecord that) {
        if (this == that) return true;
        if (that == null) return false;

        if (mSequenceIndex != that.mSequenceIndex) return false;
        if (mSequenceLength != that.mSequenceLength) return false;
        if (mSequenceName != null ? !mSequenceName.equals(that.mSequenceName) : that.mSequenceName != null)
            return false;
        // If one record has an optional attribute and the other does not, that is not considered inequality.
        
        if (mAttributes != null) {
            for (final Map.Entry<String, Object> entry: getAttributes()) {
                final Object thatAttribute = that.getAttribute(entry.getKey());
                if (thatAttribute != null && !entry.getValue().equals(thatAttribute)) {
                    return false;
                }
            }
        }

        return true;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof SAMSequenceRecord)) return false;

        final SAMSequenceRecord that = (SAMSequenceRecord) o;

        if (mSequenceIndex != that.mSequenceIndex) return false;
        if (mSequenceLength != that.mSequenceLength) return false;
        if (mAttributes != null ? !mAttributes.equals(that.mAttributes) : that.mAttributes != null) return false;
        if (mSequenceName != null ? !mSequenceName.equals(that.mSequenceName) : that.mSequenceName != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = mSequenceName != null ? mSequenceName.hashCode() : 0;
        result = 31 * result + mSequenceIndex;
        result = 31 * result + mSequenceLength;
        result = 31 * result + (mAttributes != null ? mAttributes.hashCode() : 0);
        return result;
    }
}

