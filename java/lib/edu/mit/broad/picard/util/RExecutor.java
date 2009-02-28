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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.io.IoUtil;

/**
 * Util class for executing R scripts.
 * 
 * @author Doug Voet
 */
public class RExecutor {
    private static final String R_EXE = "Rscript";
    
    /**
     * Executes the given R script that is stored in a file on the classpath. The script file
     * is read from the classpath and written to a temp file then executed by a call to Rscript.
     * Blocks until the R script is complete.
     * 
     * @param rScriptName the fully qualified name of the classpath resource of the script
     * @param arguments any arguments required by the script
     * @return the return code of the R process
     */
    public static int executeFromClasspath(String rScriptName, String... arguments) {
        File scriptFile = writeScriptFile(rScriptName);
        int returnCode = executeFromFile(scriptFile, arguments);
        scriptFile.delete();
        return returnCode;
    }

    /**
     * Executes the given R script that is stored in a file by a call to Rscript.
     * Blocks until the R script is complete.
     * 
     * @param scriptFile the file object for the script
     * @param arguments any arguments required by the script
     * @return the return code of the R process
     */
    public static int executeFromFile(File scriptFile, String... arguments) {
        String[] command = new String[arguments.length + 2];
        command[0] = R_EXE;
        command[1] = scriptFile.getAbsolutePath();
        System.arraycopy(arguments, 0, command, 2, arguments.length);
        return ProcessExecutor.execute(command);
    }

    /**
     * Writes the classpath resource named by rScriptName to the temp dir.
     */
    private static File writeScriptFile(String rScriptName) {
        InputStream scriptStream = null;
        OutputStream scriptFileStream = null;
        try {
            scriptStream = RExecutor.class.getClassLoader().getResourceAsStream(rScriptName);
            if (scriptStream == null) {
                throw new IllegalArgumentException("Script [" + rScriptName + "] not found in classpath");
            }
            File scriptFile = File.createTempFile("script", ".R");
            scriptFileStream = IoUtil.openFileForWriting(scriptFile);
            IoUtil.copyStream(scriptStream, scriptFileStream);
            return scriptFile;
        } catch (IOException e) {
            throw new PicardException("Unexpected exception creating R script file", e);
        } finally {
            if (scriptStream != null) {
                try {
                    scriptStream.close();
                } catch (IOException e) {
                }
            }
            if (scriptFileStream != null) {
                try {
                    scriptFileStream.close();
                } catch (IOException e) {
                }
            }
        }
    }
}
