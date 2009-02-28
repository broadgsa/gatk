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

import edu.mit.broad.dcp.CallStatus;

public class DistributedCallMessage
    extends DistributedMessage
{
    public DistributedCallMessage() {
    }

    public Long getCallId() {
        return mCallId;
    }

    public void setCallId(Long value) {
        mCallId = value;
    }

    public CallStatus getCallStatus() {
        return mCallStatus;
    }

    public void setCallStatus(CallStatus value) {
        mCallStatus = value;
    }

    public String getMethodName() {
        return mMethodName;
    }

    public void setMethodName(String value) {
        mMethodName = value;
    }

    public Object[] getMethodArgs() {
        return mMethodArgs;
    }

    public void setMethodArgs(Object[] value) {
        mMethodArgs = value;
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("DistributedCallMessage");
        builder.append("(");
        builder.append("" + getSenderWorkerId());
        builder.append(",");
        builder.append("" + getSenderProcessId());
        builder.append(",");
        builder.append("" + getReceiverWorkerId());
        builder.append(",");
        builder.append("" + getReceiverProcessId());
        builder.append(",");
        builder.append("" + mCallId);
        builder.append(",");
        builder.append("" + mCallStatus);
        builder.append(",");
        builder.append("" + mMethodName);
        builder.append(",");
        if (mMethodArgs == null) {
            builder.append("" + mMethodArgs);
        } else {
            builder.append("[");
            for (int i = 0; i < mMethodArgs.length; i++) {
                if (i > 0) {
                    builder.append(",");
                }
                builder.append("" + mMethodArgs[i]);
            }
            builder.append("]");
        }
        builder.append(")");
        return builder.toString();
    }

    public Long mCallId;
    public CallStatus mCallStatus;
    public String mMethodName;
    public Object[] mMethodArgs;
}
