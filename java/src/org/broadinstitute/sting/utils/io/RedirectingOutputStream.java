package org.broadinstitute.sting.utils.io;

import java.io.OutputStream;
import java.io.IOException;
/**
 * User: hanna
 * Date: Apr 30, 2009
 * Time: 5:53:32 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A stream that allows redirection to a variety of sources transparently to the
 * user of the class.
 */
@Deprecated
public class RedirectingOutputStream extends OutputStream {

    /**
     * Informs us which output stream should be used.
     */
    private OutputStreamProvider provider;

    /**
     * Build a new output stream, given the function telling us where to
     * send output.
     * @param provider Function which returns an output stream.
     */
    public RedirectingOutputStream( OutputStreamProvider provider ) {
        this.provider = provider;
    }

    /**
     * Gets the OutputStream backing this redirector now.
     * Note that the backing output stream could change at any time.
     * Use sparingly (for testing).
     */
    public OutputStream getBackingOutputStream() {
        return provider.getOutputStream();
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
     * Provides whatever output stream this data should go to at the moment.
     */
    public interface OutputStreamProvider {
        public OutputStream getOutputStream();
    }
}
