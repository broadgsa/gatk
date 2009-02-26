/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2007 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.dcp;

import java.util.*;

/**
 * Command line driver for distributed worker invocation.
 */
public class DistributedAlgorithmWorker
{
    public static void main(String[] args)
        throws Exception {
        new DistributedAlgorithmWorker().run(args);
    }

    private void run(String[] args)
        throws Exception {

        if (!parseArguments(args)) {
            System.exit(1);
        }
        System.out.println("# DistributedAlgorithmWorker");
        System.out.println("# Started at " + new Date());
        runDistributedWorker();
        System.out.println("# Ended at " + new Date());
    }

    private boolean parseArguments(String[] args) {

        int argpos = 0;
        int argsleft = 0;

        while (argpos < args.length) {
            argsleft = args.length - argpos;
            String arg = args[argpos];
            if (arg.equals("-serverHost") && argsleft > 1) {
                argpos++;
                mServerHost = args[argpos++];
            } else if (arg.equals("-serverPort") && argsleft > 1) {
                argpos++;
                mServerPort = Integer.parseInt(args[argpos++]);
            } else if (arg.equals("-workerId") && argsleft > 1) {
                argpos++;
                mWorkerId = new Integer(args[argpos++]);
            } else if (arg.equals("-processId") && argsleft > 1) {
                argpos++;
                mProcessId = new Integer(args[argpos++]);
            } else if (arg.equals("-debug")) {
                argpos++;
                mDebug = true;
                continue;
            } else if (arg.equals("-verbose")) {
                argpos++;
                mVerbose = true;
                continue;
            } else if (arg.startsWith("-")) {
                usage();
                return false;
            } else {
                break;
            }
        }

        argsleft = args.length - argpos;
        if (argsleft != 0) {
            usage();
            return false;
        }

        return true;
    }

    private void usage() {
        System.out.println("Usage: DistributedWorkerMain ...");
        System.out.println("  -serverHost <hostname>");
        System.out.println("  -serverPort <port>");
        System.out.println("  -workerId <id>");
        System.out.println("  -processId <id>");
        System.out.println("  -verbose");
        System.out.println("  -debug");
    }

    private void runDistributedWorker()
        throws Exception {

        DistributedAlgorithm algorithm = null;
        String serverAddress = getServerHost() + ":" + getServerPort();
        try {
            String url = "rmi://" + serverAddress + "/DistributedCallService";
            DistributedCallService server =
                (DistributedCallService) java.rmi.Naming.lookup(url);
            algorithm = server.getAlgorithm();
        } catch (java.rmi.ConnectException exc) {
            System.out.println("# Server " + serverAddress + " not responding.");
            return;
        }

        algorithm.setServerHost(getServerHost());
        algorithm.setServerPort(getServerPort());
        algorithm.runWorker(getWorkerId(), getProcessId());
    }

    private Integer getWorkerId() {
        return mWorkerId;
    }

    private Integer getProcessId() {
        return mProcessId;
    }

    private String getServerHost() {
        return mServerHost;
    }

    private int getServerPort() {
        return mServerPort;
    }


    private boolean mDebug = false;
    private boolean mVerbose = false;
    private String mServerHost = null;
    private int mServerPort = 0;
    private Integer mWorkerId = null;
    private Integer mProcessId = null;
}
