/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;

import edu.mit.broad.picard.PicardException;

/**
 * Utility class that will execute sub processes via Runtime.getRuntime().exec(...) and read
 * off the output from stderr and stdout of the sub process. This implementation uses a different
 * thread to read each stream: the current thread for stdout and another, internal thread for 
 * stderr. This utility is able to handle concurrent executions, spawning as many threads as
 * are required to handle the concurrent load.
 *
 * @author Doug Voet
 */
public class ProcessExecutor {
    private static Log log = Log.getInstance(ProcessExecutor.class);
    private static ExecutorService executorService = Executors.newCachedThreadPool(new ThreadFactory() {
        @Override
        public Thread newThread(Runnable r) {
            return new Thread(r, "ProcessExecutor Thread");
        }
    });
    
    /**
     * Executes the command via Runtime.getRuntime().exec() then writes stderr to log.error
     * and stdout to log.info and blocks until the command is complete.
     * 
     * @see Runtime#exec(String)
     * 
     * @param command command string
     * @return return code of command
     */
    public static int execute(String command) {
        try {
            Process process = Runtime.getRuntime().exec(command);
            return readStreamsAndWaitFor(process);
        } catch (Throwable t) {
            throw new PicardException("Unexpected exception executing [" + StringUtil.join(" ", command) + "]", t);
        }
    }

    /**
     * Executes the command via Runtime.getRuntime().exec() then writes stderr to log.error
     * and stdout to log.info and blocks until the command is complete.
     * 
     * @see Runtime#exec(String[])
     * 
     * @param commandParts command string
     * @return return code of command
     */
    public static int execute(String[] commandParts) {
        try {
            Process process = Runtime.getRuntime().exec(commandParts);
            return readStreamsAndWaitFor(process);
        } catch (Throwable t) {
            throw new PicardException("Unexpected exception executing [" + StringUtil.join(" ", commandParts) + "]", t);
        }
    }

    private static int readStreamsAndWaitFor(Process process)
            throws InterruptedException, ExecutionException {
        Future<?> stderrReader = executorService.submit(new LogErrorProcessOutputReader(process.getErrorStream()));
        new LogInfoProcessOutputReader(process.getInputStream()).run();
        // wait for stderr reader to be done
        stderrReader.get();
        return process.waitFor();
    }
    
    /**
     * Runnable that reads off the given stream and logs it somewhere.
     */
    private static abstract class ProcessOutputReader implements Runnable {
        private BufferedReader reader;
        public ProcessOutputReader(InputStream stream) {
            reader = new BufferedReader(new InputStreamReader(stream));
        }

        @Override
        public void run() {
            try {
                String line;
                while ((line = reader.readLine()) != null) {
                    log(line);
                }
            } catch (IOException e) {
                throw new PicardException("Unexpected exception reading from process stream", e);
            }
        }
        
        protected abstract void log(String message);
    }
    
    private static class LogErrorProcessOutputReader extends ProcessOutputReader {
        public LogErrorProcessOutputReader(InputStream stream) { super(stream); }
        @Override protected void log(String message) { log.error(message); }
    }

    private static class LogInfoProcessOutputReader extends ProcessOutputReader {
        public LogInfoProcessOutputReader(InputStream stream) { super(stream); }
        @Override protected void log(String message) { log.info(message); }
    }
}
