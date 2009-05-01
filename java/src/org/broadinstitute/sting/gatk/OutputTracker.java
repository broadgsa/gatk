package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.io.RedirectingOutputStream;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintWriter;
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
    protected OutputStream globalOut;
    protected OutputStream globalErr;

    /**
     * Thread-local versions of the output stream.
     */
    protected ThreadLocal<OutputStream> localOut = new ThreadLocal<OutputStream>();
    protected ThreadLocal<OutputStream> localErr = new ThreadLocal<OutputStream>();

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
            globalOut = globalErr = outputFile;
        }
        else {
            globalOut = (outFileName != null) ? prepareOutputFile( outFileName ) : System.out;
            globalErr = (errFileName != null) ? prepareOutputFile( errFileName ) : System.err;
        }
    }

    /**
     * Gets the output stream for the walker.
     * @return Output stream; should be either file-backed or System.out.
     */
    public OutputStream getGlobalOutStream() {
        return globalOut;
    }

    /**
     * Gets the error stream for the walker.
     * @return Error stream; should be either file-backed or System.err.
     */
    public OutputStream getGlobalErrStream() {
        return globalErr;
    }

    /**
     * Retrieve an output stream that will always return the most appropriate
     * writer for this thread.
     */
    public OutputStream getOutStream() {
        // Create an anonymous inner class which will just return the most
        // appropriate OutputStream from those streams stored in this class.
        return new RedirectingOutputStream(
                new RedirectingOutputStream.OutputStreamProvider() {
                    public OutputStream getOutputStream() {
                        return localOut.get() != null ? localOut.get() : globalOut;
                    }
                } );
    }

    /**
     * Retrieve the most appropriate output stream for this thread.
     * Will retrieve thread-local if available; otherwise, it'll read the
     * global stream.
     */
    public OutputStream getErrStream() {
        // Create an anonymous inner class which will just return the most
        // appropriate OutputStream from those streams stored in this class.
        return new RedirectingOutputStream(
                new RedirectingOutputStream.OutputStreamProvider() {
                    public OutputStream getOutputStream() {
                        return localErr.get() != null ? localErr.get() : globalErr;
                    }
                } );
    }

    /**
     * Set the (thread-)local version of the given output and error files.
     * @param out output stream to which to write.
     * @param err error stream to which to write.
     */
    public void setLocalStreams( OutputStream out, OutputStream err ) {
        localOut.set( out );
        localErr.set( err );
    }

    /**
     * Remove pointers to alternate, local output streams.
     */
    public void removeLocalStreams() {
        localOut.remove();
        localErr.remove();
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
