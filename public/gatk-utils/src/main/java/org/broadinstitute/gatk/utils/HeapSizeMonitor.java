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

package org.broadinstitute.gatk.utils;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;

/**
 * Monitor the current heap size, allowing the application to programmatically
 * access the data.
 *
 * @author mhanna
 * @version 0.1
 */
public class HeapSizeMonitor {
    private final int monitorFrequencyMillis;
    private final MonitorRunnable monitorRunnable;

    private Thread monitorThread;

    public HeapSizeMonitor() {
        this(1000);
    }

    public HeapSizeMonitor(final int monitorFrequencyMillis) {
        this.monitorFrequencyMillis = monitorFrequencyMillis;
        this.monitorRunnable = new MonitorRunnable();
    }

    public long getMaxMemoryUsed() {
        return monitorRunnable.getMaxMemoryUsed();
    }

    public void start() {
        monitorThread = new Thread(monitorRunnable);
        monitorThread.start();
    }

    public void stop() {
        monitorRunnable.stop = true;
        try {
            monitorThread.join();
        }
        catch(InterruptedException ex) {
            throw new ReviewedGATKException("Unable to connect to monitor thread");
        }
        monitorThread = null;        
    }

    private class MonitorRunnable implements Runnable {
        private MemoryMXBean monitor;

        private long maxMemoryUsed;
        private boolean stop;

        public MonitorRunnable() {
            monitor = ManagementFactory.getMemoryMXBean();   
        }

        public void reset() {
            maxMemoryUsed = 0L;
            stop = false;
        }

        public long getMaxMemoryUsed() {
            return maxMemoryUsed;
        }

        public void run() {
            while(!stop) {
                System.gc();
                maxMemoryUsed = Math.max(monitor.getHeapMemoryUsage().getUsed(),maxMemoryUsed);
                try {
                    Thread.sleep(monitorFrequencyMillis);
                }
                catch(InterruptedException ex) {
                    throw new ReviewedGATKException("Unable to continue monitoring heap consumption",ex);
                }
            }
        }
    }
}
