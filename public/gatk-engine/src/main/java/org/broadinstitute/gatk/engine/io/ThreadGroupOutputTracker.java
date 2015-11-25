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

package org.broadinstitute.gatk.engine.io;

import org.broadinstitute.gatk.engine.executive.OutputMergeTask;
import org.broadinstitute.gatk.engine.io.storage.Storage;
import org.broadinstitute.gatk.engine.io.storage.StorageFactory;
import org.broadinstitute.gatk.engine.io.stubs.Stub;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * An output tracker that can either track its output per-thread or directly.
 *
 * This output tracker doesn't use thread local values, but rather looks up the
 * storage map via the thread's group.  This is necessary in the case where
 * there's a master thread that creates the output map, and spawns subthreads
 * that actually do work.  As long as those subthreads are spawned in the
 * thread group of the master thread, this tracker will properly find the
 * storage map associated with the master thread in the group, and return
 * the map to all subthreads.
 *
 * @author mhanna, depristo
 * @version 0.2
 */
public class ThreadGroupOutputTracker extends OutputTracker {
    /**
     * A map from thread ID of the master thread to the storage map from
     * Stub to Storage objects
     */
    private Map<ThreadGroup, Map<Stub, Storage>> threadsToStorage = new HashMap<ThreadGroup, Map<Stub, Storage>>();

    /**
     * A total hack.  If bypass = true, bypass thread local storage and write directly
     * to the target file.  Used to handle output during initialize() and onTraversalDone().
     */
    private boolean bypass = false;
    public void bypassThreadLocalStorage(boolean bypass) {
        this.bypass = bypass;
    }

    /**
     * Initialize the storage map for this thread.
     *
     * Checks if there's a thread local binding for this thread, and if
     * not initializes the map for it.  This map is then
     * populated with stub -> storage bindings according to the
     * superclasses' outputs map.
     *
     * Must be called within the master thread to create a map associated with
     * the master thread ID.
     */
    public synchronized void initializeStorage() {
        final ThreadGroup group = Thread.currentThread().getThreadGroup();
        Map<Stub,Storage> threadLocalOutputStreams = threadsToStorage.get(group);

        if( threadLocalOutputStreams == null ) {
            threadLocalOutputStreams = new HashMap<Stub,Storage>();
            threadsToStorage.put( group, threadLocalOutputStreams );
        }

        for ( final Stub stub : outputs.keySet() ) {
            final Storage target = StorageFactory.createStorage(stub, createTempFile(stub));
            threadLocalOutputStreams.put(stub, target);
        }
    }

    @Override
    public <T> T getStorage( final Stub<T> stub ) {
        Storage target;

        if (bypass) {
            target = outputs.get(stub);
            if( target == null ) {
                target = StorageFactory.createStorage(stub);
                outputs.put(stub, target);
            }
        }
        else {
            final Map<Stub,Storage> threadLocalOutputStreams = findStorage(Thread.currentThread());
            target = threadLocalOutputStreams.get(stub);

            // make sure something hasn't gone wrong, and we somehow find a map that doesn't include our stub
            if ( target == null )
                throw new ReviewedGATKException("target isn't supposed to be null for " + Thread.currentThread()
                        + " id " + Thread.currentThread().getId() + " map is " + threadLocalOutputStreams);
        }

        return (T)target;
    }


    private synchronized Map<Stub,Storage> findStorage(final Thread thread) {
        final Map<Stub, Storage> map = threadsToStorage.get(thread.getThreadGroup());

        if ( map != null ) {
            return map;
        } else {
            // something is terribly wrong, we have a storage lookup for a thread that doesn't have
            // any map data associated with it!
            throw new ReviewedGATKException("Couldn't find storage map associated with thread " + thread + " in group " + thread.getThreadGroup());
        }
    }

    /**
     * Close down any existing temporary files which have been opened.
     */
    public synchronized OutputMergeTask closeStorage() {
        final Map<Stub,Storage> threadLocalOutputStreams = findStorage(Thread.currentThread());

        if( threadLocalOutputStreams == null || threadLocalOutputStreams.isEmpty() )
            return null;

        final OutputMergeTask outputMergeTask = new OutputMergeTask();
        for( Map.Entry<Stub,Storage> entry: threadLocalOutputStreams.entrySet() ) {
            final Stub stub = entry.getKey();
            final Storage storageEntry = entry.getValue();

            storageEntry.close();
            outputMergeTask.addMergeOperation(getTargetStream(stub), storageEntry);
        }

//        logger.info("Closing " + Thread.currentThread().getId() + " => " + threadLocalOutputStreams);
        threadLocalOutputStreams.clear();

        return outputMergeTask;
    }

    /**
     * Creates a temporary file for a stub of the given type.
     * @param stub Stub for which to create a temporary file.
     * @param <T> Type of the stub to accept.
     * @return A temp file, or throw an exception if the temp file cannot be created.
     */
    private <T> File createTempFile( Stub<T> stub ) {
        try {
            return File.createTempFile( stub.getClass().getName(), null );
        } catch( IOException ex ) {
            throw new UserException.BadTmpDir("Unable to create temporary file for stub: " + stub.getClass().getName() );
        }
    }
}
