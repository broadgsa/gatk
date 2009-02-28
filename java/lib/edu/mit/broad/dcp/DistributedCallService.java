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

public interface DistributedCallService
    extends java.rmi.Remote
{
    public DistributedAlgorithm getAlgorithm()
        throws java.rmi.RemoteException;
    public long writeMessage(DistributedCallMessage message)
        throws java.rmi.RemoteException;
    public DistributedCallMessage acceptMessage(int workerId, int processId)
        throws java.rmi.RemoteException;
    public void completeMessage(int workerId, int processId, long callId)
        throws java.rmi.RemoteException;
}
