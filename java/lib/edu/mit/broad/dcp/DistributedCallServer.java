/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.dcp;


import edu.mit.broad.dcp.message.*;

import java.rmi.server.UnicastRemoteObject;
import java.util.*;

public class DistributedCallServer
    extends UnicastRemoteObject
    implements DistributedCallService
{
    public DistributedCallServer()
        throws java.rmi.RemoteException {
    }

    public void setAlgorithm(DistributedAlgorithm algorithm) {
        mAlgorithm = algorithm;
    }

    public DistributedAlgorithm getAlgorithm() {
        return mAlgorithm;
    }

    public long writeMessage(DistributedCallMessage message) {
        message.setCallStatus(CallStatus.PENDING);
        message.setCallId(generateCallId());
        if (message.getReceiverWorkerId().equals(0)) {
            synchronized (mMessageQueue) {
                mMessageQueue.addLast(message);
            }
        } else {
            synchronized (mMessageQueue) {
                mMessageQueue.addFirst(message);
            }
        }            
        return message.getCallId();
    }

    public DistributedCallMessage acceptMessage(int workerId, int processId) {
        if (workerId <= 0) {
            throw new IllegalArgumentException("Invalid worker ID: " + workerId);
        }
        if (processId <= 0) {
            throw new IllegalArgumentException("Invalid process ID: " + processId);
        }
        synchronized (mMessageQueue) {
            Iterator<DistributedCallMessage> iterator = mMessageQueue.iterator();
            while (iterator.hasNext()) {
                DistributedCallMessage message = iterator.next();
                if (message.getCallStatus() != CallStatus.PENDING) {
                    continue;
                }
                int receiverId = message.getReceiverWorkerId();
                if (receiverId == workerId ||
                    (receiverId == 0 && workerId > 1)) {
                    message.setCallStatus(CallStatus.PROCESSING);
                    message.setReceiverWorkerId(workerId);
                    message.setReceiverProcessId(processId);
                    return message;
                }
            }
        }

        return null;
    }

    public void completeMessage(int workerId, int processId, long callId) {
        if (workerId <= 0) {
            throw new IllegalArgumentException("Invalid worker ID: " + workerId);
        }
        if (processId <= 0) {
            throw new IllegalArgumentException("Invalid process ID: " + processId);
        }
        if (callId <= 0) {
            throw new IllegalArgumentException("Invalid call ID: " + callId);
        }
        synchronized (mMessageQueue) {
            Iterator<DistributedCallMessage> iterator = mMessageQueue.iterator();
            while (iterator.hasNext()) {
                DistributedCallMessage message = iterator.next();
                if (message.getCallId().longValue() == callId) {
                    if (message.getCallStatus() != CallStatus.PROCESSING) {
                        throw new IllegalStateException("Call #" + callId + " not in state PROCESSING");
                    }
                    if (!message.getReceiverWorkerId().equals(workerId)) {
                        throw new IllegalStateException("Call #" + callId + " assigned to worker " + message.getReceiverWorkerId() + " not worker " + workerId);
                    }
                    if (!message.getReceiverProcessId().equals(processId)) {
                        throw new IllegalStateException("Call #" + callId + " assigned to process " + message.getReceiverProcessId() + " not process " + processId);
                    }
                    iterator.remove();
                    return;
                }
            }
        }

        throw new IllegalArgumentException("Unrecognized call ID " + callId);
    }

    public boolean isQueueEmpty() {
        synchronized (mMessageQueue) {
            return mMessageQueue.isEmpty();
        }
    }

    public void stop() {
        try {
            UnicastRemoteObject.unexportObject(this, false);
        } catch (java.rmi.NoSuchObjectException exc) {
            throw new RuntimeException("Exception unexporting object: " + exc.getMessage(),
                                       exc);
        }
    }

    private synchronized long generateCallId() {
        return ++mCallIdGenerator;
    }

    private long mCallIdGenerator = 0;
    private DistributedAlgorithm mAlgorithm = null;
    private LinkedList<DistributedCallMessage> mMessageQueue =
        new LinkedList<DistributedCallMessage>();
}
