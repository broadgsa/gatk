/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.io;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import edu.mit.broad.picard.PicardException;

/**
 * A class for utility methods that wrap or aggregate functionality in Java IO.
 *
 * @author Tim Fennell
 */
public class IoUtil {
    /**
     * Checks that a file is non-null, exists, is not a directory and is readable.  If any
     * condition is false then a runtime exception is thrown.
     *
     * @param file the file to check for readability
     */
    public static void assertFileIsReadable(File file) {
        if (file == null) {
			throw new IllegalArgumentException("Cannot check readability of null file.");
		} else if (!file.exists()) {
            throw new PicardException("Cannot read non-existent file: " + file.getAbsolutePath());
        }
        else if (file.isDirectory()) {
            throw new PicardException("Cannot read file because it is a directory: " + file.getAbsolutePath());
        }
        else if (!file.canRead()) {
            throw new PicardException("File exists but is not readable: " + file.getAbsolutePath());
        }
    }

    /**
     * Checks that a file is non-null, and is either extent and writable, or non-existent but
     * that the parent directory exists and is writable. If any
     * condition is false then a runtime exception is thrown.
     *
     * @param file the file to check for writability
     */
    public static void assertFileIsWritable(File file) {
        if (file == null) {
			throw new IllegalArgumentException("Cannot check readability of null file.");
		} else if (!file.exists()) {
            // If the file doesn't exist, check that it's parent directory does and is writable
            File parent = file.getAbsoluteFile().getParentFile();
            if (!parent.exists()) {
                throw new PicardException("Cannot write file: " + file.getAbsolutePath() + ". " +
                        "Neither file nor parent directory exist.");
            }
            else if (!parent.isDirectory()) {
                throw new PicardException("Cannot write file: " + file.getAbsolutePath() + ". " +
                        "File does not exist and parent is not a directory.");
            }
            else if (!parent.canWrite()) {
                throw new PicardException("Cannot write file: " + file.getAbsolutePath() + ". " +
                        "File does not exist and parent directory is not writable..");
            }
        }
        else if (file.isDirectory()) {
            throw new PicardException("Cannot write file because it is a directory: " + file.getAbsolutePath());
        }
        else if (!file.canWrite()) {
            throw new PicardException("File exists but is not writable: " + file.getAbsolutePath());
        }
    }

    /**
     * Checks that a directory is non-null, extent, writable and a directory 
     * otherwise a runtime exception is thrown.
     *
     * @param dir the dir to check for writability
     */
    public static void assertDirectoryIsWritable(File dir) {
        if (dir == null) {
            throw new IllegalArgumentException("Cannot check readability of null file.");
        } 
        else if (!dir.exists()) {
            throw new PicardException("Directory does not exist: " + dir.getAbsolutePath());
        }
        else if (!dir.isDirectory()) {
            throw new PicardException("Cannot write to directory because it is not a directory: " + dir.getAbsolutePath());
        }
        else if (!dir.canWrite()) {
            throw new PicardException("Directory exists but is not writable: " + dir.getAbsolutePath());
        }
    }

    /**
     * Opens a file for reading, decompressing it if necessary
     *
     * @param file  The file to open
     * @return the input stream to read from
     */
    public static InputStream openFileForReading(File file) {

        try {
            if (file.getName().endsWith(".gz") ||
                file.getName().endsWith(".bfq") ||
                file.getName().endsWith(".map")) {
                return new GZIPInputStream(new FileInputStream(file));
            }
            //TODO: Other compression formats
            else {
                return new FileInputStream(file);
            }
        }
        catch (IOException ioe) {
            throw new PicardException("File not found: " + file.getName(), ioe);
        }

    }

    /**
     * Opens a file for writing, overwriting the file if it already exists
     *
     * @param file  the file to write to
     * @return the output stream to write to
     */
    public static OutputStream openFileForWriting(File file) {
        return openFileForWriting(file, false);
    }

    /**
     * Opens a file for writing
     *
     * @param file  the file to write to
     * @param append    whether to append to the file if it already exists (we overwrite it if false)
     * @return the output stream to write to
     */
    public static OutputStream openFileForWriting(File file, boolean append) {

        try {
            if (file.getName().endsWith(".gz") ||
                file.getName().endsWith(".bfq") ||
                file.getName().endsWith(".map")) {
                return new GZIPOutputStream(new FileOutputStream(file, append));
            }
            //TODO: Other compression formats
            else {
                return new FileOutputStream(file, append);
            }
        }
        catch (IOException ioe) {
            throw new PicardException("Error opening file for writing: " + file.getName(), ioe);
        }
    }
    
    /**
     * Utility method to copy the contents of input to output. The caller is responsible for
     * opening and closing both streams.
     * 
     * @param input contents to be copied
     * @param output destination
     */
    public static void copyStream(InputStream input, OutputStream output) {
        try {
            byte[] buffer = new byte[1024];
            int bytesRead = 0;
            while((bytesRead = input.read(buffer)) > 0) {
                output.write(buffer, 0, bytesRead);
            }
        } catch (IOException e) {
            throw new PicardException("Exception copying stream", e);
        }
    }

}
