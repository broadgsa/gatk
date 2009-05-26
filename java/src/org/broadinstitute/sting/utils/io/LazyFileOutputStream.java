package org.broadinstitute.sting.utils.io;

import org.broadinstitute.sting.utils.StingException;

import java.io.OutputStream;
import java.io.IOException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
/**
 * User: hanna
 * Date: May 26, 2009
 * Time: 3:51:49 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * An output stream that only initializes itself the first time its used.
 * Needs a callback that can create an output stream.
 */
public class LazyFileOutputStream extends OutputStream {
    /**
     * Generates output files on demand.
     */
    private final FileFactory factory;

    private File targetFile = null;

    /**
     * The target for any writes performed by the output stream.
     */
    private FileOutputStream targetOutputStream = null;

    /**
     * Create a new LazyOutputStream, indicating how to create a new stream.
     * @param factory Creator of the output stream, when necessary.
     */
    public LazyFileOutputStream( FileFactory factory ) {
        this.factory = factory;
    }

    /**
     * Indicates whether the LazyOutputStream had to get off its butt and create
     * a new output stream.
     * @return
     */
    public boolean isCreated() {
        return targetOutputStream != null;
    }

    /**
     * Public method to return the lazily created file.
     * @return Stream created by the lazy loader.
     * @throw StingException if no stream was created.
     */
    public File getBackingFile() {
        if( targetFile == null )
            throw new StingException("No lazy-loaded stream was created.");
        return targetFile;
    }

    @Override
    public void close() throws IOException {
        getBackingOutputStream().close();
    }

    @Override
    public void flush() throws IOException {
        getBackingOutputStream().flush();
    }

    @Override
    public void write(byte[] b) throws IOException {
        getBackingOutputStream().write(b);
    }

    @Override
    public void write(byte[] b, int off, int len) throws IOException {
        getBackingOutputStream().write(b,off,len);
    }

    @Override
    public void write(int b) throws IOException {
        getBackingOutputStream().write(b);
    }

    /**
     * Lazy loader for the output stream.
     */
    protected OutputStream getBackingOutputStream() {
        if( targetOutputStream == null ) {
            try {
                targetFile = factory.create();
                targetOutputStream = new FileOutputStream( targetFile );
            }
            catch( IOException ex ) {
                throw new StingException("Unable to open new temp file", ex );
            }
        }
        return targetOutputStream;
    }        

    /**
     * Teaches the LazyOutputStream how to create a new outputstream when necessary.
     */
    public interface FileFactory {
        public File create() throws IOException;
    }
}
