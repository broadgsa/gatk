package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.io.LazyFileOutputStream;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.nio.channels.FileChannel;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;

/**
 * User: hanna
 * Date: Apr 30, 2009
 * Time: 4:04:38 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Hold pointers to the output and error streams, and state to indicate whether
 * a write is complete.  Not generally thread-safe.  Calls to isComplete()/complete()
 * can be made at any time from any thread, but complete() should be called on the
 * thread which is doing the writing. 
 */
public class OutputMerger {
    /**
     * Is writing to these streams complete?
     */
    private boolean complete = false;

    /**
     * The printstreams which should be written to.
     */
    private LazyFileOutputStream out = null;
    private LazyFileOutputStream err = null;

    public OutputMerger() {
        out = new LazyFileOutputStream( outStreamFactory );
        err = new LazyFileOutputStream( errStreamFactory );
    }

    /**
     * Creates a generic output stream for temporary data.
     */
    private static class TempFileFactory implements LazyFileOutputStream.FileFactory {
        private final String prefix;
        public TempFileFactory( String prefix ) { this.prefix = prefix; }
        public File create() throws IOException { return File.createTempFile(prefix,null); }
    }

    /**
     * Creates a stream for temporary out (not err) data.
     */
    private static final TempFileFactory outStreamFactory = new TempFileFactory("gatkout_");

    /**
     * Creates a stream for temporary err data.
     */
    private static final TempFileFactory errStreamFactory = new TempFileFactory("gatkerr_");

    /**
     * Waits for any the given OutputMerger to be ready for merging.
     */
    public synchronized void waitForOutputComplete() {
        try {
            wait();
        }
        catch( InterruptedException ex ) {
            throw new StingException("Interrupted while waiting for more output to be finalized.",ex);            
        }
    }

    /**
     * Is this output complete?
     * @return True if output complete.  False otherwise.
     */
    public synchronized boolean isComplete() {
        return complete;
    }

    /**
     * Indicate that no more data will be written to these output streams.
     */
    public synchronized void complete() {
        if( isComplete() )
            throw new IllegalStateException("Tried to complete an output merge twice.");

        try {
            if( out.isCreated() ) {
                out.flush();
                out.close();
            }
            if( err.isCreated() ) {
                err.flush();
                err.close();
            }
        }
        catch( IOException ex ) {
            throw new StingException( "Unable to close sharding output files", ex );
        }

        this.complete = true;

        // Notify waiting listeners that this shard is complete and ready for merging.
        notifyAll();
    }

    public void mergeInto( OutputStream globalOut, OutputStream globalErr ) {
        synchronized(this) {
            if( !isComplete() )
                throw new StingException("Tried to merge incomplete stream into output");
        }

        if( out.isCreated() ) transferContentsToTarget( out.getBackingFile(), globalOut );
        if( err.isCreated() ) transferContentsToTarget( err.getBackingFile(), globalErr );
    }

    /**
     * Copy the contents of the given file into the specified output stream.
     * @param source Source of data to copy.
     * @param target Target for copied data.
     */
    private void transferContentsToTarget( File source, OutputStream target ) {
        FileInputStream sourceStream = null;
        try {
            sourceStream = new FileInputStream( source );
            FileChannel sourceChannel = sourceStream.getChannel();

            WritableByteChannel targetChannel = Channels.newChannel( target );
            sourceChannel.transferTo( 0, sourceChannel.size(), targetChannel );

            sourceStream.close();
            source.delete();            
        }
        catch( FileNotFoundException ex ) {
            throw new StingException("Unable to open input stream for file: " + source,ex);
        }
        catch( IOException ex ) {
            throw new StingException("Unable to transfer contents of file:" + source,ex);
        }
    }

    /**
     * Return the stream where output is sent.
     * @return output stream.
     */
    public OutputStream getOutStream() {
        return out;
    }

    /**
     * Gets the stream where error info is sent.
     * @return error stream.
     */
    public OutputStream getErrStream() {
        return err;
    }
}
