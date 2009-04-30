package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.StingException;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.FileNotFoundException;
/**
 * User: hanna
 * Date: Apr 30, 2009
 * Time: 9:40:09 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Manages the output and err streams that are created specifically for walker
 * output.
 */

public class OutputTracker {
    /**
     * The streams to which walker users should be writing directly.
     */
    private PrintStream out;
    private PrintStream err;

    /**
     * Cache the file streams so that reassemblers can use nio to do fast file transfers.
     */
    private FileOutputStream outFileStream;
    private FileOutputStream errFileStream;


    /**
     * Create an object to manage output given filenames for the output and error files.
     * If no files are specified, returns null.
     * @param outFileName Name of the output file.
     * @param errFileName Name of the error file.
     */
    public OutputTracker( String outFileName, String errFileName ) {
        // If the two output streams match and are non-null, initialize them identically.
        // Otherwise, initialize them separately.
        if( outFileName != null && outFileName.equals(errFileName) ) {
            FileOutputStream outputFile = prepareOutputFile( outFileName );
            outFileStream = errFileStream = outputFile;
            out = err = new PrintStream( outputFile );
        }
        else {
            if( outFileName != null ) {
                outFileStream = prepareOutputFile( outFileName );
                out = new PrintStream( outFileStream );
            }
            else
                out = System.out;

            if( errFileName != null ) {
                errFileStream = prepareOutputFile( errFileName );
                err = new PrintStream( errFileStream );
            }
            else
                err = System.err;
        }
    }

    /**
     * Gets the output stream for the walker.
     * @return Output stream; should be either file-backed or System.out.
     */
    public PrintStream getOutStream() {
        return out;
    }

    /**
     * Gets the error stream for the walker.
     * @return Error stream; should be either file-backed or System.err.
     */
    public PrintStream getErrStream() {
        return err;
    }

    /**
     * Gets the filestream associated with normal output.
     * @return FileStream associated with the output; null if not backed by a file.
     */
    public FileOutputStream getOutFile() {
        return outFileStream;
    }

    /**
     * Gets the filestream associated with error output.
     * @return stream associated with error output; null if not backed by a file.
     */
    public FileOutputStream getErrFile() {
        return errFileStream;
    }

    /**
     * Given a (non-null) filename, open a file for output.
     * @param fileName Filename for the stream.  Should not be null.
     * @return An open output stream associated with the given file.
     */
    private FileOutputStream prepareOutputFile( String fileName ) {
        FileOutputStream outputFile = null;

        try {
            outputFile = new FileOutputStream( fileName );
        }
        catch( FileNotFoundException ex ) {
            throw new StingException("Unable to open output file " + fileName, ex);
        }

        return outputFile;
    }
}
