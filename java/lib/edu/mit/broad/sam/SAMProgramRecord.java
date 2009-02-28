/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class SAMProgramRecord {
    public static final String PROGRAM_GROUP_ID_TAG = "ID";
    private static final String PROGRAM_VERSION_TAG = "VN";
    private static final String COMMAND_LINE_TAG = "CL";
    private final String mProgramGroupId;
    private final Map<String, String> mAttributes = new HashMap<String, String>();

    public SAMProgramRecord(final String programGroupId) {
        this.mProgramGroupId = programGroupId;
    }

    public String getProgramGroupId() {
        return mProgramGroupId;
    }

    public String getAttribute(final String key) {
        return mAttributes.get(key);
    }

    public void setAttribute(final String key, final String value) {
        mAttributes.put(key, value);
    }

    public Set<Map.Entry<String, String>> getAttributes() {
        return mAttributes.entrySet();
    }

    public String getProgramVersion() {
        return getAttribute(PROGRAM_VERSION_TAG);
    }

    public void setProgramVersion(final String version) {
        setAttribute(PROGRAM_VERSION_TAG, version);
    }

    public String getCommandLine() {
        return getAttribute(COMMAND_LINE_TAG);
    }

    public void setCommandLine(final String commandLine) {
        setAttribute(COMMAND_LINE_TAG, commandLine);
    }

    /**
     * @return true if this == that except for the program group ID, which is arbitrary
     */
    public boolean equivalent(final SAMProgramRecord that) {
        return mAttributes.equals(that.mAttributes);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SAMProgramRecord that = (SAMProgramRecord) o;

        if (mAttributes != null ? !mAttributes.equals(that.mAttributes) : that.mAttributes != null) return false;
        if (mProgramGroupId != null ? !mProgramGroupId.equals(that.mProgramGroupId) : that.mProgramGroupId != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = mProgramGroupId != null ? mProgramGroupId.hashCode() : 0;
        result = 31 * result + (mAttributes != null ? mAttributes.hashCode() : 0);
        return result;
    }
}
