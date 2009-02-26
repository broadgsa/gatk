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

import edu.mit.broad.dcp.message.*;

import java.io.*;
import java.util.*;
import java.lang.reflect.Method;
import java.net.InetAddress;
import java.net.ServerSocket;
import java.rmi.registry.*;

/**
 * Experimental.
 */
public abstract class DistributedAlgorithm
    implements Serializable
{
    public static final Integer ANY = 0;
    public static final Integer MASTER = 1;

    public DistributedAlgorithm() {
    }

    public String getServerHost() {
        return mServerHost;
    }

    public void setServerHost(String value) {
        mServerHost = value;
    }

    public int getServerPort() {
        return mServerPort;
    }

    public void setServerPort(int value) {
        mServerPort = value;
    }

    public String getAlgorithmName() {
        if (mAlgorithmName != null) {
            return mAlgorithmName;
        } else {
            return getClassName();
        }
    }

    public void setAlgorithmName(String value) {
        mAlgorithmName = value;
    }

    public int getMaximumWorkerCount() {
        return mMaximumWorkerCount;
    }

    public void setMaximumWorkerCount(int value) {
        mMaximumWorkerCount = value;
    }

    /**
     * Name of LSF queue to use for workers.
     */
    public String getLsfQueue() {
        return mLsfQueue;
    }

    public void setLsfQueue(String value) {
        mLsfQueue = value;
    }

    /**
     * Directory to hold lsf log files.
     */
    public String getLsfLogDirectory() {
        return mLsfLogDirectory;
    }

    public void setLsfLogDirectory(String value) {
        mLsfLogDirectory = value;
    }

    public boolean getEnableGcLogging() {
        return mEnableGcLogging;
    }

    public void setEnableGcLogging(boolean value) {
        mEnableGcLogging = value;
    }

    public Integer getWorkerId() {
        return mWorkerId;
    }

    public Integer getProcessId() {
        return mProcessId;
    }

    protected void init()
        throws Exception {
    }

    protected abstract void start()
        throws Exception;

    public void run()
        throws Exception {

        if (mIsRunning) {
            throw new IllegalStateException("Algorithm is already running");
        }

        mIsRunning = true;
        mWorkerId = MASTER;
        mProcessId = MASTER;

        try {
            startDistributedServer();
            init();
            startWorkerThread();
            startWorkers();
            start();
            waitForCompletion();
        } finally {
            // TBD: More cleanup (shutdown threads, etc.)
            stopDistributedServer();
            mIsRunning = false;
        }
    }

    void runWorker(int workerId, int processId)
        throws Exception {

        if (mIsRunning) {
            throw new IllegalStateException("Algorithm is already running");
        }

        mIsRunning = true;
        mWorkerId = workerId;
        mProcessId = processId;

        try {
            if (openDistributedServer() == null) {
                report("Server " + mServerHost + ":" + mServerPort + " not responding");
                return;
            }
            init();
            startWorkerThread();
            mWorkerThread.join();
        } finally {
            closeDistributedServer();
            mIsRunning = false;
        }
    }

    private void startWorkers() {
        int workerCount = getMaximumWorkerCount();
        if (workerCount <= 0) {
            // Use single process execution for testing/debugging.
            new InProcessWorker().start();
            return;
        }
        if (workerCount > 1000) {
            throw new RuntimeException("Excessive worker count: " + workerCount);
        }
        for (int i = 0; i < workerCount; i++) {
            Integer workerId = (MASTER + i + 1);
            Integer processId = workerId;  // for now
            startWorker(workerId, processId);
        }
    }

    private void startDistributedServer() {
        try {
            // Create a server socket to allocate a unique port.
            // There is a window of vulnerability where the port
            // can get reused, but in practice this works ok.
            String serverHost = getCurrentHost();
            ServerSocket socket = new ServerSocket(0);
            int serverPort = socket.getLocalPort();
            socket.close();
            Registry registry = LocateRegistry.createRegistry(serverPort);
            DistributedCallServer server = new DistributedCallServer();
            server.setAlgorithm(this);
            registry.bind("DistributedCallService", server);
            mServerHost = serverHost;
            mServerPort = serverPort;
            mDistributedCallServer = server;
            mDistributedCallService = server;
        } catch (Exception exc) {
            throw wrapException(exc);
        }
    }

    private void stopDistributedServer() {
        if (mDistributedCallServer != null) {
            try {
                Registry registry = LocateRegistry.getRegistry(mServerPort);
                registry.unbind("DistributedCallService");
                mDistributedCallServer.stop();
            } catch (Exception exc) {
                throw wrapException(exc);
            }
        }
        mDistributedCallService = null;
        mDistributedCallServer = null;
    }

    private DistributedCallService openDistributedServer() {
        mDistributedCallService = null;
        try {
            String url = "rmi://" + getServerHost() + ":" + getServerPort() + "/DistributedCallService";
            DistributedCallService server =
                (DistributedCallService) java.rmi.Naming.lookup(url);
            mDistributedCallService = server;
        } catch (java.rmi.NotBoundException exc) {
            // Server has exited
        } catch (Exception exc) {
            throw wrapException(exc);
        }
        return mDistributedCallService;
    }

    private void closeDistributedServer() {
        mDistributedCallService = null;
    }

    private void startWorker(Integer workerId, Integer processId) {

        String logFile = "worker_" + processId + "_%J.bsub";
        if (mLsfLogDirectory != null) {
            logFile = mLsfLogDirectory + "/" + logFile;
        }

        List<String> command = new ArrayList<String>();
        command.add("bsub");
        command.add("-o");
        command.add(logFile);
        if (mLsfQueue != null) {
            command.add("-q");
            command.add(mLsfQueue);
        }
        command.add("runDistributedWorker");
        command.add("-serverHost");
        command.add(getServerHost());
        command.add("-serverPort");
        command.add(Integer.toString(getServerPort()));
        command.add("-workerId");
        command.add(Integer.toString(workerId));
        command.add("-processId");
        command.add(Integer.toString(processId));

        // Pass our -Xmx setting along to all workers.
        Map<String, String> environment =
            new LinkedHashMap<String, String>(System.getenv());
        long maxMemory = Runtime.getRuntime().maxMemory();
        long maxKbytes = maxMemory / 1024;
        String memJavaOpt = "-Xmx" + maxKbytes + "K";

        // Enable GC logging if requested
        String gcJavaOpt = null;
        if (mEnableGcLogging) {
            String gcLogFile = "worker_" + processId + ".gc.log";
            if (mLsfLogDirectory != null) {
                gcLogFile = mLsfLogDirectory + "/" + gcLogFile;
            }
            gcJavaOpt = "-Xloggc:" + gcLogFile;
        }

        String javaOpts = environment.get("JAVAOPTS");
        if (javaOpts == null) {
            javaOpts = memJavaOpt;
            if (gcJavaOpt != null) {
                javaOpts = javaOpts + " " + gcJavaOpt;
            }
            environment.put("JAVAOPTS", javaOpts);
        }

        // Log output ourselves (rather than waiting for bsub).
        String workerLogFile = "worker_" + processId + ".log";
        if (mLsfLogDirectory != null) {
            workerLogFile = mLsfLogDirectory + "/" + workerLogFile;
        }
        environment.put("DA_LOG_FILE", workerLogFile);

        CommandRunner runner = new CommandRunner();
        Writer output = new LsfOutputFilter();
        runner.setStandardOutputDestination(output);
        runner.setStandardErrorDestination(output);
        String[] commandArray = command.toArray(new String[command.size()]);
        String[] environmentArray = createEnvironmentArray(environment);
        int status = runner.runCommand(commandArray, environmentArray, null);
        if (status != 0) {
            throw new RuntimeException("Error starting worker: " + status);
        }
    }

    private String[] createEnvironmentArray(Map<String, String> map) {
        if (map == null) {
            return null;
        }
        int index = 0;
        String[] array = new String[map.size()];
        for (Map.Entry<String, String> entry : map.entrySet()) {
            array[index++] = entry.getKey() + "=" + entry.getValue();
        }
        return array;
    }

    private String getCurrentHost() {
        try {
            return InetAddress.getLocalHost().getCanonicalHostName();
        } catch (Exception exc) {
            throw wrapException(exc);
        }
    }

    private void waitForCompletion() {
        DistributedCallServer server = mDistributedCallServer;
        while (true) {
            if (server.isQueueEmpty()) {
                break;
            }
            try {
                Thread.sleep(1000);
            } catch (InterruptedException exc) {
                // ignore
            }
        }
    }

    protected void callDistributed(String methodName, Object... methodArgs) {
        callDistributed(null, methodName, methodArgs);
    }

    protected void callDistributed(Integer workerId, String methodName, Object... methodArgs) {
        if (workerId == null) {
            workerId = ANY;
        }
        try {
            DistributedCallMessage message = new DistributedCallMessage();
            message.setSenderWorkerId(getWorkerId());
            message.setSenderProcessId(getProcessId());
            message.setReceiverWorkerId(workerId);
            message.setMethodName(methodName);
            message.setMethodArgs(methodArgs);
            mDistributedCallService.writeMessage(message);
        } catch (Throwable exc) {
            throw wrapException(exc);
        }
    }

    private void callMethod(String methodName, Object[] methodArgs) {
        try {
            Object target = this;
            Class targetClass = target.getClass();
            Method targetMethod = findMethod(targetClass, methodName);
            if (targetMethod == null) {
                throw new RuntimeException("Cannot find target method: " + methodName);
            }
            targetMethod.invoke(target, methodArgs);
        } catch (Throwable exc) {
            throw wrapException(exc);
        }
    }

    private Method findMethod(Class clazz, String methodName) throws Exception {
        Method result = null;
        Method[] methods = clazz.getDeclaredMethods();
        for (int i = 0; i < methods.length; i++) {
            if (methods[i].getName().equals(methodName)) {
                if (result != null) {
                    throw new RuntimeException("Duplicate method name: " + methodName);
                }
                result = methods[i];
            }
        }
        return result;
    }

    private RuntimeException wrapException(Throwable exception) {
        if (exception instanceof RuntimeException) {
            return (RuntimeException) exception;
        } else {
            return new RuntimeException(exception.getMessage(), exception);
        }
    }

    private void startWorkerThread() {
        if (mWorkerThread != null) {
            throw new IllegalStateException("WorkerThread is running");
        }
        mWorkerThread = new WorkerThread();
        mWorkerThread.start();
    }

    private void stopWorkerThread() {
        if (mWorkerThread == null) {
            throw new IllegalStateException("WorkerThread is running");
        }
        mWorkerThread.stopThread();
    }

    private class WorkerThread extends Thread {

        WorkerThread() {
            setDaemon(true);
        }

        public void run() {
            try {
                DistributedCallService service = mDistributedCallService;
                while (true) {
                    if (isInterrupted()) {
                        System.out.println("#DBG: Worker isInterrupted");
                        throw new InterruptedException();
                    }
                    DistributedCallMessage message =
                        service.acceptMessage(getWorkerId(), getProcessId());
                    if (message == null) {
                        Thread.sleep(1000);
                    } else {
                        processMessage(message);
                    }
                }
            } catch (InterruptedException exc) {
                // Interruption terminates this thread.
                // System.out.println("#DBG: Worker caught InterruptedException");
            } catch (Throwable exc) {
                if (isDisconnectException(exc)) {
                    report("Server disconnected");
                } else {
                    reportError("Exception in WorkerThread: " + exc.getMessage(), exc);
                    System.exit(1);
                }
            }
            report("WorkerThread terminated");
        }

        void stopThread() {
            // System.out.println("#DBG: About to interrupt worker...");
            interrupt();
            // System.out.println("#DBG: Joining worker...");
            try {
                join();
            } catch (InterruptedException exc) {
                // ignore
            }
        }

        private boolean isDisconnectException(Throwable exc) {
            if (exc instanceof java.rmi.ConnectException) {
                return true;
            } else if (exc instanceof java.rmi.NoSuchObjectException) {
                return true;
            } else if (exc instanceof java.rmi.UnmarshalException &&
                       exc.getCause() != null &&
                       exc.getCause() instanceof EOFException) {
                return true;
            } else {
                return false;
            }
        }
    }

    private void processMessage(DistributedCallMessage message) {
        try {
            Integer workerId = message.getReceiverWorkerId();
            if (workerId == null || !workerId.equals(getWorkerId())) {
                reportError("Invalid worker ID in message: " + message);
                return;
            }
            callMethod(message.getMethodName(), message.getMethodArgs());
        } catch (Throwable exc) {
            reportError("Exception running message: " + message, exc);
        } finally {
            completeMessage(message);
        }
    }

    private void completeMessage(DistributedCallMessage message) {
        try {
            DistributedCallService service = mDistributedCallService;
            service.completeMessage(getWorkerId(), getProcessId(), message.getCallId());
        } catch (Throwable exc) {
            reportError("Exception completing message: " + message, exc);
        }
    }

    protected void report(String message) {
        String identity =
            getAlgorithmName() + " " +
            getWorkerId() + "/" + getProcessId();
        System.out.println("# " + identity + " : " + message);
    }

    protected void reportError(String message) {
        reportError(message, null);
    }

    protected void reportError(String message, Throwable exception) {
        String identity =
            getAlgorithmName() + " " +
            getWorkerId() + "/" + getProcessId();
        System.out.println("Error" +
                           " [" + identity + "]" +
                           ": " + message);
        if (exception != null) {
            System.out.println(" with exception: " + exception.getMessage());
            exception.printStackTrace(System.out);
        }
    }

    private String getClassName() {
        String name = getClass().getName();
        return name.substring(name.lastIndexOf('.')+1);
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("DistributedAlgorithm");
        builder.append("(");
        builder.append("" + getAlgorithmName());
        builder.append(",");
        builder.append("" + getWorkerId());
        builder.append(",");
        builder.append("" + getProcessId());
        builder.append(",");
        builder.append("" + getMaximumWorkerCount());
        builder.append(",");
        builder.append("" + getLsfQueue());
        builder.append(",");
        builder.append("" + mIsRunning);
        builder.append(")");
        return builder.toString();
    }

    // This class is used only during in-process execution/testing/debugging.
    private class InProcessWorker extends Thread {

        InProcessWorker() {
            setDaemon(true);
        }

        public void run() {
            report("InProcessWorker starting");
            try {
                String serverAddress = getServerHost() + ":" + getServerPort();
                String url = "rmi://" + serverAddress + "/DistributedCallService";
                DistributedCallService server =
                    (DistributedCallService) java.rmi.Naming.lookup(url);
                DistributedAlgorithm algorithm = server.getAlgorithm();
                algorithm.setServerHost(getServerHost());
                algorithm.setServerPort(getServerPort());
                algorithm.runWorker(2, 1);
            } catch (Throwable exc) {
                reportError("Exception in InProcessWorker: " + exc.getMessage(), exc);
                System.exit(1);
            }
            report("InProcessWorker terminated");
        }
    }

    private static class LsfOutputFilter
        extends FilterWriter {

        LsfOutputFilter() {
            super(new PrintWriter(System.out, true));
        }

        public void write(int ch)
            throws IOException {
            if (mAtLineStart) {
                out.write("# ");
                mAtLineStart = false;
            }
            out.write(ch);
            mAtLineStart = (ch == '\n');
        }

        public void write(String s, int off, int len)
            throws IOException {
            write(s.toCharArray(), off, len);
        }

        public void write(char[] a, int off, int len)
            throws IOException {
            for (int i = 0; i < len; i++) {
                write(a[off+i]);
            }
        }

        private boolean mAtLineStart = true;
    }


    private transient int mMaximumWorkerCount = 0;
    private transient String mLsfQueue = null;
    private transient String mLsfLogDirectory = null;
    private transient boolean mEnableGcLogging = false;
    private transient boolean mIsRunning = false;
    private transient int mWorkerId = 0;
    private transient int mProcessId = 0;
    private transient WorkerThread mWorkerThread = null;
    private transient String mAlgorithmName = null;
    private transient String mServerHost = null;
    private transient int mServerPort = 0;
    private transient DistributedCallService mDistributedCallService = null;
    private transient DistributedCallServer mDistributedCallServer = null;
}
