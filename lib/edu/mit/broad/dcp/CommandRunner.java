/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2006 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.dcp;

import java.io.*;


/**
 * Utility class to run system commands synchronously and return the output.
 * 
 * The interface supports the typical case where you want to return a modest
 * amount of information from the command's standard output or standard error
 * as a string.  The caller can override this behavior, however, and provide
 * alternative output destinations if necessary.
 * 
 * If setMergeOutput() is true, then this class will attempt to interleave
 * the standard output and standard error streams of the command into one
 * stream (standard output).  This may not produce exactly the same results
 * as having the operating system interleave the output, but works well for
 * simple executables that do not heavily intermix stdout and stderr.
 * 
 * A typical invocation is:
 * <pre>
 *  CommandRunner runner = new CommandRunner();
 *  int status = runner.runCommand("ls");
 *  if (status == 0) {
 *      System.out.print(runner.getStandardOutput());
 *  }
 * </pre>
 * 
 * @author Bob Handsaker
 */
public class CommandRunner {

    private boolean mMergeOutput = false;
    private Writer mStandardOutputDestination = null;
    private Writer mStandardErrorDestination = null;
    private String mStandardOutputString = null;
    private String mStandardErrorString = null;


    /**
     * Default constructor.
     */
    public CommandRunner() {
    }

    /**
     * Get the standard output from the last command as a string.
     * 
     * If no command has been run or an explicit output destination
     * was set, then this method returns null.
     */
    public String getStandardOutputString() {
        return mStandardOutputString;
    }

    /**
     * Get the standard error from the last command as a string.
     * 
     * If no command has been run or an explicit output destination
     * was set, then this method returns null.
     */
    public String getStandardErrorString() {
        return mStandardErrorString;
    }

    /**
     * If true, the command's standard error stream will be interleaved
     * with the command's standard output stream.  The standard error
     * stream destination will not be used.
     */
    public boolean getMergeOutput() {
        return mMergeOutput;
    }

    /**
     * If true, the command's standard error stream will be interleaved
     * with the command's standard output stream.
     */
    public void setMergeOutput(boolean value) {
        mMergeOutput = value;
    }

    /**
     * The destination for the command's standard output stream.
     * If null, the standard output will be captured in a string.
     */
    public Writer getStandardOutputDestination() {
        return mStandardOutputDestination;
    }

    /**
     * The destination for the command's standard output stream.
     * If set to null, the standard output will be captured in a string.
     */
    public void setStandardOutputDestination(Writer writer) {
        mStandardOutputDestination = writer;
    }

    /**
     * The destination for the command's standard error stream.
     * If null, the standard error will be captured in a string.
     */
    public Writer getStandardErrorDestination() {
        return mStandardErrorDestination;
    }

    /**
     * The destination for the command's standard error stream.
     * If set to null, the standard error will be captured in a string.
     */
    public void setStandardErrorDestination(Writer writer) {
        mStandardErrorDestination = writer;
    }

    /**
     * Run a command string as a system command.
     * 
     * Returns the exit status of the command.
     * 
     * When this method is called, the standard output string
     * and standard error string are updated if no alternative output
     * destinations have been set.
     * 
     * This method throws a RuntimeException if running the command fails
     * (for example, if there are not enough system resources to spawn
     * the process).
     * 
     * @param commmand The command string to run.
     * @return Command exit status.
     * @throws RuntimeException If command execution fails.
     */
    public int runCommand(String command)
        throws RuntimeException {
        return runCommand(command.split(" "), null, null);
    }

    /**
     * Run a command string as a system command.
     * 
     * Returns the exit status of the command.
     * 
     * When this method is called, the standard output string
     * and standard error string are updated if no alternative output
     * destinations have been set.
     * 
     * This method throws a RuntimeException if running the command fails
     * (for example, if there are not enough system resources to spawn
     * the process).
     * 
     * @param commmand The command string to run.
     * @param environment The command environment (or null to inherit).
     * @param workingDirectory The working directory (or null to inherit).
     * @return Command exit status.
     * @throws RuntimeException If command execution fails.
     */
    public int runCommand(String command, String[] environment, File workingDirectory)
        throws RuntimeException {
        return runCommand(command.split(" "), environment, workingDirectory);
    }

    /**
     * Run a command string as a system command.
     * 
     * Returns the exit status of the command.
     * 
     * When this method is called, the standard output string
     * and standard error string are updated if no alternative output
     * destinations have been set.
     * 
     * This method throws a RuntimeException if running the command fails
     * (for example, if there are not enough system resources to spawn
     * the process).
     * 
     * @param commmand The command to run (as a array of arguments).
     * @param environment The command environment (or null to inherit).
     * @param workingDirectory The working directory (or null to inherit).
     * @return Command exit status.
     * @throws RuntimeException If command execution fails.
     */
    public int runCommand(String[] command, String[] environment, File workingDirectory)
        throws RuntimeException {

        Writer stdout = mStandardOutputDestination;
        Writer stderr = mStandardErrorDestination;
        if (stdout == null) {
            stdout = new StringWriter();
        }
        if (mMergeOutput) {
            stderr = stdout;
        } else if (stderr == null) {
            stderr = new StringWriter();
        }

        mStandardOutputString = null;
        mStandardErrorString = null;

        int commandStatus = 0;
        try {
            Process process =
                Runtime.getRuntime().exec(command, environment, workingDirectory);
            StreamHandler stdoutHandler =
                new StreamHandler(process.getInputStream(), stdout);
            StreamHandler stderrHandler =
                new StreamHandler(process.getErrorStream(), stderr);
            
            commandStatus = process.waitFor();

            // Wait for the streams to drain.
            stdoutHandler.join();
            stderrHandler.join();
        } catch (Exception exc) {
            throw new RuntimeException("Command execution failed: " +
                                       exc.getMessage(),
                                       exc);
        }

        if (mStandardOutputDestination == null) {
            mStandardOutputString = stdout.toString();
        }
        if (mStandardErrorDestination == null && !mMergeOutput) {
            mStandardErrorString =  stderr.toString();
        }

        return commandStatus;
    }


    /**
     * Internal class to asynchronously read from the standard output
     * and standard error streams of the command being executed.
     * 
     * If you do not handle command output asynchronously, then execution
     * of a command may block in some environments if the program produces
     * too much output.  In this case, the call to run the process will
     * never complete.
     */
    private static class StreamHandler extends Thread {

        /**
         * Constructor.
         * Create an instance of this class, which is an asynchronous
         * thread that will consume input from the given input stream
         * and send the output to the given output destination.
         *
         * @param input The input stream to read.
         * @param output The output destination.
         */
        StreamHandler(InputStream input, Writer output) {
            m_input = input;
            m_output = output;
            start();
        }


        /**
         * Standard thread run method.
         * Pipe input from the input source to the output destination
         * until there is no more input left.
         *
         * If an IOException occurs, the thread will make sure all
         * available output has been flushed to the destination and
         * then terminate.  The IOException is not propagated.
         */
        public void run() {

            char[] buffer = new char[4096];
            Reader reader =
                new InputStreamReader(new BufferedInputStream(m_input));
            Writer writer = m_output;

            try {
                while (true) {
                    int count = reader.read(buffer);
                    if (count <= 0) {
                        break;
                    }
                    if (writer != null) {
                        synchronized (writer) {
                            writer.write(buffer, 0, count);
                        }
                    }
                }
            } catch (IOException ignore) {
                // Ignore IO exceptions
            } finally {
                try {
                    reader.close();
                } catch (Exception ignore) {
                }
                try {
                    m_output.flush();
                } catch (Exception ignore) {
                }
            }
        }

        private InputStream m_input;
        private Writer m_output;
    }
}
