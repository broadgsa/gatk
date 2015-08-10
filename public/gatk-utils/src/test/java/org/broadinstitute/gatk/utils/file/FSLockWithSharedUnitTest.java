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

package org.broadinstitute.gatk.utils.file;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.io.File;

public class FSLockWithSharedUnitTest extends BaseTest {

    private static final int MAX_EXPECTED_LOCK_ACQUISITION_TIME = FSLockWithShared.DEFAULT_LOCK_ACQUISITION_TIMEOUT_IN_MILLISECONDS +
                                                                  FSLockWithShared.THREAD_TERMINATION_TIMEOUT_IN_MILLISECONDS;

    /**
     * Test to ensure that we're never spending more than the maximum configured amount of time in lock acquisition calls.
     */
    @Test( timeOut = MAX_EXPECTED_LOCK_ACQUISITION_TIME + 10 * 1000 )
    public void testLockAcquisitionTimeout() {
        final File lockFile = createTempFile("FSLockWithSharedUnitTest", ".lock");
        final FSLockWithShared lock = new FSLockWithShared(lockFile);
        boolean lockAcquisitionSucceeded = false;

        try {
            lockAcquisitionSucceeded = lock.sharedLock();
        }
        catch ( UserException e ) {
            logger.info("Caught UserException from lock acquisition call: lock acquisition must have timed out. Message: " + e.getMessage());
        }
        finally {
            if ( lockAcquisitionSucceeded ) {
                lock.unlock();
            }
        }
    }
}
