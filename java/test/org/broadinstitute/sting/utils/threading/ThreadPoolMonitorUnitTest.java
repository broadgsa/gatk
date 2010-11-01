package org.broadinstitute.sting.utils.threading;

import org.testng.annotations.Test;
import org.broadinstitute.sting.BaseTest;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors; 
/**
 * User: hanna
 * Date: Apr 29, 2009
 * Time: 4:30:55 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Tests for the thread pool monitor class.
 */

public class ThreadPoolMonitorUnitTest extends BaseTest {
    private ExecutorService threadPool = Executors.newFixedThreadPool(1);

    /**
     * Test to make sure the thread pool wait works properly. 
     */
    @Test(timeOut=2000)
    public void testThreadPoolMonitor() {
        ThreadPoolMonitor monitor = new ThreadPoolMonitor();
        synchronized(monitor) {
            threadPool.submit(monitor);
            monitor.watch();
        }
    }
}
