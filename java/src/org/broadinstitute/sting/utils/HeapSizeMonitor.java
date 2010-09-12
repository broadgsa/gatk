package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

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
            throw new ReviewedStingException("Unable to connect to monitor thread");
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
                    throw new ReviewedStingException("Unable to continue monitoring heap consumption",ex);
                }
            }
        }
    }
}
