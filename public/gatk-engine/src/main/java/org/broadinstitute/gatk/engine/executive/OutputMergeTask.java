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

package org.broadinstitute.gatk.engine.executive;

import org.broadinstitute.gatk.engine.io.storage.Storage;

import java.util.ArrayList;
import java.util.Collection;

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
public class OutputMergeTask {
    /**
     * The output streams which should be written to.
     */
    private final Collection<MergeOperation<?>> mergeOperations = new ArrayList<MergeOperation<?>>();

    /**
     * Add a new merge operation to this merge task.
     * @param targetStream Target for stream output.
     * @param temporaryStorage Temporary storage.
     * @param <StreamType> Type of the output stream.
     */
    public <StreamType> void addMergeOperation( StreamType targetStream, Storage<StreamType> temporaryStorage ) {
        mergeOperations.add( new MergeOperation<StreamType>(targetStream,temporaryStorage) );
    }

    /**
     * Merge data from output streams into target storage.
     */
    public synchronized void merge() {
        for( MergeOperation mergeOperation: mergeOperations )
            mergeOperation.temporaryStorage.mergeInto(mergeOperation.targetStream);
    }

    /**
     * Represents a single file needed to be merged.
     * @param <StreamType> Type of the file to be merged.
     */
    private class MergeOperation<StreamType> {
        /**
         * Destination for the temporary file's output.
         */
        public final StreamType targetStream;

        /**
         * Temporary storage location for the file.
         */
        public final Storage<StreamType> temporaryStorage;

        /**
         * Create a new merge file object with the given output stream and storage placeholder.
         * @param targetStream Target for temporary data.
         * @param temporaryStorage The temporary data itself.
         */
        public MergeOperation( StreamType targetStream, Storage<StreamType> temporaryStorage ) {
            this.targetStream = targetStream;
            this.temporaryStorage = temporaryStorage;
        }
    }
}

