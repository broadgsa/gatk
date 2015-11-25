/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.jna.drmaa.v1_0;

import com.sun.jna.Pointer;
import org.ggf.drmaa.*;

import java.util.*;

/**
 * JNA mapping from Java to C DRMAA binding.
 */
public class JnaJobTemplate implements JobTemplate {
    private final JnaSession session;
    private final Pointer jt;

    public JnaJobTemplate(JnaSession session, Pointer jt) {
        this.session = session;
        this.jt = jt;
    }

    public Pointer getPointer() {
        return jt;
    }

    @Override
    public void setRemoteCommand(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_REMOTE_COMMAND, s);
    }

    @Override
    public String getRemoteCommand() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_REMOTE_COMMAND);
    }

    @SuppressWarnings("unchecked")
    @Override
    public void setArgs(List list) throws DrmaaException {
        JnaSession.setVectorAttribute(jt, LibDrmaa.DRMAA_V_ARGV, list);
    }

    @Override
    public List getArgs() throws DrmaaException {
        return JnaSession.getVectorAttribute(jt, LibDrmaa.DRMAA_V_ARGV);
    }

    @Override
    public void setJobSubmissionState(int state) throws DrmaaException {
        String stateString;
        if (state == JobTemplate.HOLD_STATE)
            stateString = LibDrmaa.DRMAA_SUBMISSION_STATE_HOLD;
        else if (state == JobTemplate.ACTIVE_STATE)
            stateString = LibDrmaa.DRMAA_SUBMISSION_STATE_ACTIVE;
        else
            throw new InvalidAttributeValueException("jobSubmissionState attribute is invalid");
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_JS_STATE, stateString);
    }

    @Override
    public int getJobSubmissionState() throws DrmaaException {
        int state;
        String stateString = JnaSession.getAttribute(jt, LibDrmaa.DRMAA_JS_STATE);
        if (LibDrmaa.DRMAA_SUBMISSION_STATE_HOLD.equals(stateString))
            state = JobTemplate.HOLD_STATE;
        else if (LibDrmaa.DRMAA_SUBMISSION_STATE_ACTIVE.equals(stateString))
            state = JobTemplate.ACTIVE_STATE;
        else
            throw new InvalidAttributeValueException("jobSubmissionState attribute is invalid");
        return state;
    }

    @SuppressWarnings("unchecked")
    @Override
    public void setJobEnvironment(Map env) throws DrmaaException {
        JnaSession.setVectorAttribute(jt, LibDrmaa.DRMAA_V_ENV, JnaSession.mapToCollection(env));
    }

    @SuppressWarnings("unchecked")
    @Override
    public Map getJobEnvironment() throws DrmaaException {
        return JnaSession.collectionToMap(JnaSession.getVectorAttribute(jt, LibDrmaa.DRMAA_V_ENV));
    }

    @Override
    public void setWorkingDirectory(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_WD, s);
    }

    @Override
    public String getWorkingDirectory() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_WD);
    }

    @Override
    public void setJobCategory(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_JOB_CATEGORY, s);
    }

    @Override
    public String getJobCategory() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_JOB_CATEGORY);
    }

    @Override
    public void setNativeSpecification(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_NATIVE_SPECIFICATION, s);
    }

    @Override
    public String getNativeSpecification() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_NATIVE_SPECIFICATION);
    }

    @SuppressWarnings("unchecked")
    @Override
    public void setEmail(Set set) throws DrmaaException {
        JnaSession.setVectorAttribute(jt, LibDrmaa.DRMAA_V_EMAIL, set);
    }

    @SuppressWarnings("unchecked")
    @Override
    public Set getEmail() throws DrmaaException {
        return new LinkedHashSet<String>(JnaSession.getVectorAttribute(jt, LibDrmaa.DRMAA_V_EMAIL));
    }

    @Override
    public void setBlockEmail(boolean b) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_BLOCK_EMAIL, b ? "1" : "0");
    }

    @Override
    public boolean getBlockEmail() throws DrmaaException {
        return "1".equals(JnaSession.getAttribute(jt, LibDrmaa.DRMAA_BLOCK_EMAIL));
    }

    @Override
    public void setStartTime(PartialTimestamp partialTimestamp) throws DrmaaException {
        JnaSession.setPartialTime(jt, LibDrmaa.DRMAA_START_TIME, partialTimestamp);
    }

    @Override
    public PartialTimestamp getStartTime() throws DrmaaException {
        return JnaSession.getPartialTime(jt, LibDrmaa.DRMAA_START_TIME);
    }

    @Override
    public void setJobName(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_JOB_NAME, s);
    }

    @Override
    public String getJobName() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_JOB_NAME);
    }

    @Override
    public void setInputPath(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_INPUT_PATH, s);
    }

    @Override
    public String getInputPath() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_INPUT_PATH);
    }

    @Override
    public void setOutputPath(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_OUTPUT_PATH, s);
    }

    @Override
    public String getOutputPath() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_OUTPUT_PATH);
    }

    @Override
    public void setErrorPath(String s) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_ERROR_PATH, s);
    }

    @Override
    public String getErrorPath() throws DrmaaException {
        return JnaSession.getAttribute(jt, LibDrmaa.DRMAA_ERROR_PATH);
    }

    @Override
    public void setJoinFiles(boolean b) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_JOIN_FILES, b ? "y" : "n");
    }

    @Override
    public boolean getJoinFiles() throws DrmaaException {
        return "y".equals(JnaSession.getAttribute(jt, LibDrmaa.DRMAA_JOIN_FILES));
    }

    @Override
    public void setTransferFiles(FileTransferMode fileTransferMode) throws DrmaaException {
        StringBuilder buf = new StringBuilder();

        if (fileTransferMode.getInputStream())
            buf.append('i');

        if (fileTransferMode.getOutputStream())
            buf.append('o');

        if (fileTransferMode.getErrorStream())
            buf.append('e');

        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_TRANSFER_FILES, buf.toString());
    }

    @Override
    public FileTransferMode getTransferFiles() throws DrmaaException {
        String mode = JnaSession.getAttribute(jt, LibDrmaa.DRMAA_TRANSFER_FILES);

        if (mode == null)
            return null;

        FileTransferMode fileTransferMode = new FileTransferMode();
        fileTransferMode.setInputStream(mode.indexOf('i') >= 0);
        fileTransferMode.setOutputStream(mode.indexOf('o') >= 0);
        fileTransferMode.setErrorStream(mode.indexOf('e') >= 0);
        return fileTransferMode;
    }

    @Override
    public void setDeadlineTime(PartialTimestamp partialTimestamp) throws DrmaaException {
        JnaSession.setPartialTime(jt, LibDrmaa.DRMAA_DEADLINE_TIME, partialTimestamp);
    }

    @Override
    public PartialTimestamp getDeadlineTime() throws DrmaaException {
        return JnaSession.getPartialTime(jt, LibDrmaa.DRMAA_DEADLINE_TIME);
    }

    @Override
    public void setHardWallclockTimeLimit(long l) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_WCT_HLIMIT, JnaSession.formatLimit(l));
    }

    @Override
    public long getHardWallclockTimeLimit() throws DrmaaException {
        return JnaSession.parseLimit(JnaSession.getAttribute(jt, LibDrmaa.DRMAA_WCT_HLIMIT));
    }

    @Override
    public void setSoftWallclockTimeLimit(long l) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_WCT_SLIMIT, JnaSession.formatLimit(l));
    }

    @Override
    public long getSoftWallclockTimeLimit() throws DrmaaException {
        return JnaSession.parseLimit(JnaSession.getAttribute(jt, LibDrmaa.DRMAA_WCT_SLIMIT));
    }

    @Override
    public void setHardRunDurationLimit(long l) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_DURATION_HLIMIT, JnaSession.formatLimit(l));
    }

    @Override
    public long getHardRunDurationLimit() throws DrmaaException {
        return JnaSession.parseLimit(JnaSession.getAttribute(jt, LibDrmaa.DRMAA_DURATION_HLIMIT));
    }

    @Override
    public void setSoftRunDurationLimit(long l) throws DrmaaException {
        JnaSession.setAttribute(jt, LibDrmaa.DRMAA_DURATION_SLIMIT, JnaSession.formatLimit(l));
    }

    @Override
    public long getSoftRunDurationLimit() throws DrmaaException {
        return JnaSession.parseLimit(JnaSession.getAttribute(jt, LibDrmaa.DRMAA_DURATION_SLIMIT));
    }

    @Override
    public Set getAttributeNames() throws DrmaaException {
        return JnaSession.getAttrNames();
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof JnaJobTemplate))
            return false;
        JnaJobTemplate other = (JnaJobTemplate) obj;
        return this.jt.equals(other.jt) && this.session.equals(other.session);
    }

    @Override
    public int hashCode() {
        return jt.hashCode();
    }
}
