/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.io.stubs;

import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.engine.io.OutputTracker;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

/**
 * A stub for routing and management of anything backed by an OutputStream.
 *
 * @author mhanna
 * @version 0.1
 */
public class OutputStreamStub extends OutputStream implements Stub<OutputStream> {
    /**
     * The file that this stub should write to.  Should be passed along to
     * whatever happens to create storage for this stub.  Might be null, if
     * this stub connects directly to an existing stream.
     */
    private final File targetFile;

    /**
     * The stream that this stub should write to.  Should be passed along to
     * whatever happens to create storage for this stub.  Might be null, if
     * this stub connects directly to an existing stream.
     */
    private final OutputStream targetStream;
    
    /**
     * Connects this stub with an external stream capable of serving the
     * requests of the consumer of this stub.
     */
    private OutputTracker outputTracker = null;

    /**
     * Specify that this target output stream should write to the given file.
     * @param targetFile Target file to which to write.  Should not be null.
     */
    public OutputStreamStub( File targetFile ) {
        this.targetFile = targetFile;
        this.targetStream = null;
    }

    /**
     * Specify that this target output stream should write to the given stream.
     * @param targetStream Target stream to which to write.  Should not be null.
     */
    public OutputStreamStub( OutputStream targetStream ) {
        this.targetFile = null;
        this.targetStream = targetStream;
    }


    /**
     * Return the target file to which this data should be written.
     * @return Target file.  No sanity checking will have been performed by the file object.
     */
    public File getOutputFile() {
        return targetFile;
    }

    /**
     * Return the target stream to which this data should be written.
     * @return Target stream.  No sanity checking will have been performed by the file object.
     */
    public OutputStream getOutputStream() {
        return targetStream;
    }

    /**
     * Registers the given streamConnector with this stub.
     * @param outputTracker The connector used to provide an appropriate stream.
     */
    public void register( OutputTracker outputTracker ) {
        this.outputTracker = outputTracker;
    }

    @Override
    public void processArguments( final GATKArgumentCollection argumentCollection ) {}

    /**
     * @{inheritDoc}
     */
    public void flush() throws IOException {
        outputTracker.getStorage(this).flush();
    }

    /**
     * @{inheritDoc}
     */
    public void close() throws IOException {
        outputTracker.getStorage(this).close();
    }

    /**
     * @{inheritDoc}
     */
    public void write( byte[] b ) throws IOException {
        outputTracker.getStorage(this).write(b);
    }

    /**
     * @{inheritDoc}
     */
    public void write( byte[] b, int off, int len ) throws IOException {
        outputTracker.getStorage(this).write(b, off, len);
    }

    /**
     * @{inheritDoc}
     */
    public void write( int b ) throws IOException {
        outputTracker.getStorage(this).write(b);
    }
}
