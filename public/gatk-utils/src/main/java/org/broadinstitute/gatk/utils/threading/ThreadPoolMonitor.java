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

package org.broadinstitute.gatk.utils.threading;

import org.apache.log4j.Logger;
/**
 * User: hanna
 * Date: Apr 29, 2009
 * Time: 4:27:58 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Waits for a signal to come through that the thread pool has run
 * a given task and therefore has a free slot.
 *
 * Make sure, that, when using, the submit and the run are both
 * protected by the same synchronized(monitor) lock.  See the test
 * case for an example.
 */
public class ThreadPoolMonitor implements Runnable {
    /**
     * Logger for reporting interruptions, etc.
     */
    private static Logger logger = Logger.getLogger(ThreadPoolMonitor.class);

    /**
     * Watch the monitor
     */
    public synchronized void watch() {
        try {
            wait();
        }
        catch( InterruptedException ex ) {
            logger.error("ThreadPoolMonitor interrupted:" + ex.getStackTrace());
            throw new RuntimeException("ThreadPoolMonitor interrupted", ex);
        }
    }

    /**
     * Instruct the monitor that the thread pool has run for the class.
     * Only the thread pool should execute this method.
     */
    public synchronized void run() {
        notify();
    }
}

