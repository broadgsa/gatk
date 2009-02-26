/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2007 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.dcp.message;


public class DistributedMessage
{
    public DistributedMessage() {
    }

    public Integer getSenderWorkerId() {
        return mSenderWorkerId;
    }

    public void setSenderWorkerId(Integer value) {
        mSenderWorkerId = value;
    }

    public Integer getSenderProcessId() {
        return mSenderProcessId;
    }

    public void setSenderProcessId(Integer value) {
        mSenderProcessId = value;
    }

    public Integer getReceiverWorkerId() {
        return mReceiverWorkerId;
    }

    public void setReceiverWorkerId(Integer value) {
        mReceiverWorkerId = value;
    }

    public Integer getReceiverProcessId() {
        return mReceiverProcessId;
    }

    public void setReceiverProcessId(Integer value) {
        mReceiverProcessId = value;
    }

    public Integer mSenderWorkerId;
    public Integer mSenderProcessId;
    public Integer mReceiverWorkerId;
    public Integer mReceiverProcessId;
}
