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
 * Header information about a read group.
 */
public class SAMReadGroupRecord
{
    private String mReadGroupId = null;
    private final Map<String, Object> mAttributes = new HashMap<String, Object>();
    public static final String READ_GROUP_ID_TAG = "ID";
    public static final String READ_GROUP_SAMPLE_TAG = "SM";
    public static final String PREDICTED_MEDIAN_INSERT_SIZE_TAG = "PI";
    public static final String DATE_RUN_PRODUCED_TAG = "DT";

    public SAMReadGroupRecord(final String id) {
        mReadGroupId = id;
    }

    public String getReadGroupId() {
        return mReadGroupId;
    }

    public String getSample() {
        return (String) getAttribute("SM");
    }

    public void setSample(final String value) {
        setAttribute("SM", value);
    }

    public String getLibrary() {
        return (String) getAttribute("LB");
    }

    public void setLibrary(final String value) {
        setAttribute("LB", value);
    }

    public Object getAttribute(final String key) {
        return mAttributes.get(key);
    }

    public void setAttribute(final String key, final Object value) {
        mAttributes.put(key, value);
    }

    public Set<Map.Entry<String, Object>> getAttributes() {
        return mAttributes.entrySet();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SAMReadGroupRecord that = (SAMReadGroupRecord) o;

        if (mAttributes != null ? !mAttributes.equals(that.mAttributes) : that.mAttributes != null) return false;
        if (mReadGroupId != null ? !mReadGroupId.equals(that.mReadGroupId) : that.mReadGroupId != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = mReadGroupId != null ? mReadGroupId.hashCode() : 0;
        result = 31 * result + (mAttributes != null ? mAttributes.hashCode() : 0);
        return result;
    }
}

