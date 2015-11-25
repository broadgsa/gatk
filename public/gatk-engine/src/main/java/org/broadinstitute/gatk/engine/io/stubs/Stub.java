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
import java.io.OutputStream;

/**
 * A stub used for managing IO. Acts as a proxy for IO streams
 * not yet created or streams that need significant external
 * management.
 *
 * @author mhanna
 * @version 0.1
 */
public interface Stub<StreamType> {
    /**
     * Provides a facility to register this stream with the given
     * StreamConnector.  The stub should route each output method
     * to the specified connector.
     * @param outputTracker The connector used to provide an appropriate stream.
     */
    public void register( OutputTracker outputTracker );

    /**
     * Provides a mechanism for uniformly processing command-line arguments
     * that are important for file processing.  For example, this method
     * might pass on the compression value specified by the user to
     * a SAMFileWriter
     * @param argumentCollection The arguments to be processed
     */
    public void processArguments( final GATKArgumentCollection argumentCollection );

    /**
     * Returns the OutputStream represented by this stub or null if not available.
     */
    public OutputStream getOutputStream();

    /**
     * Returns the File represented by this stub or null if not available.
     */
    public File getOutputFile();
}
