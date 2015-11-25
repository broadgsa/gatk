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

package org.broadinstitute.gatk.utils.jna.lsf.v7_0_6;

import com.sun.jna.*;
import com.sun.jna.ptr.IntByReference;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.jna.lsf.v7_0_6.LibBat.*;

import java.io.File;

/**
 * Really unit tests, but these tests will only run on systems with LSF set up.
 */
public class LibBatQueueTest extends BaseTest {
    @BeforeClass(enabled=false)
    public void initLibBat() {
        Assert.assertFalse(LibBat.lsb_init("LibBatQueueTest") < 0, LibBat.lsb_sperror("lsb_init() failed"));
    }

    @Test(enabled=false)
    public void testClusterName() {
        String clusterName = LibLsf.ls_getclustername();
        System.out.println("Cluster name: " + clusterName);
        Assert.assertNotNull(clusterName);
    }

    @Test(enabled=false)
    public void testReadConfEnv() {
        LibLsf.config_param[] configParams = (LibLsf.config_param[]) new LibLsf.config_param().toArray(4);

        configParams[0].paramName = "LSF_UNIT_FOR_LIMITS";
        configParams[1].paramName = "LSF_CONFDIR";
        configParams[2].paramName = "MADE_UP_PARAMETER";

        Structure.autoWrite(configParams);

        if (LibLsf.ls_readconfenv(configParams[0], null) != 0) {
            Assert.fail(LibLsf.ls_sysmsg());
        }

        Structure.autoRead(configParams);

        System.out.println("LSF_UNIT_FOR_LIMITS: " + configParams[0].paramValue);
        Assert.assertNotNull(configParams[1].paramValue);
        Assert.assertNull(configParams[2].paramValue);
        Assert.assertNull(configParams[3].paramName);
        Assert.assertNull(configParams[3].paramValue);
    }

    @Test(enabled=false)
    public void testReadQueueLimits() {
        String queue = "hour";
        StringArray queues = new StringArray(new String[] {queue});
        IntByReference numQueues = new IntByReference(1);
        queueInfoEnt queueInfo = LibBat.lsb_queueinfo(queues, numQueues, null, null, 0);

        Assert.assertEquals(numQueues.getValue(), 1);
        Assert.assertNotNull(queueInfo);
        Assert.assertEquals(queueInfo.queue, queue);

        int runLimit = queueInfo.rLimits[LibLsf.LSF_RLIMIT_RUN];
        Assert.assertTrue(runLimit > 0, "LSF run limit is not greater than zero: " + runLimit);
    }

    @Test(enabled=false)
    public void testSubmitEcho() throws Exception {
        if ( ! queueTestRunModeIsSet ) {
            throw new SkipException("Skipping testSubmitEcho because we are in queue test dry run mode");
        }

        String queue = "hour";
        File outFile = tryCreateNetworkTempFile("LibBatQueueTest.out");

        submit req = new submit();

        for (int i = 0; i < LibLsf.LSF_RLIM_NLIMITS; i++)
            req.rLimits[i] = LibLsf.DEFAULT_RLIMIT;

        req.projectName = "LibBatQueueTest";
        req.options |= LibBat.SUB_PROJECT_NAME;

        req.queue = queue;
        req.options |= LibBat.SUB_QUEUE;

        req.outFile = outFile.getPath();
        req.options |= LibBat.SUB_OUT_FILE;

        req.userPriority = 100;
        req.options2 |= LibBat.SUB2_JOB_PRIORITY;

        req.command = "echo \"Hello world.\"";

        String[] argv = {"", "-a", "tv"};
        int setOptionResult = LibBat.setOption_(argv.length, new StringArray(argv), "a:", req, ~0, ~0, ~0, null);
        Assert.assertTrue(setOptionResult != -1, "setOption_ returned -1");

        submitReply reply = new submitReply();
        long jobId = LibBat.lsb_submit(req, reply);

        Assert.assertFalse(jobId < 0, LibBat.lsb_sperror("Error dispatching"));

        System.out.println("Waiting for job to run: " + jobId);
        int jobStatus = LibBat.JOB_STAT_PEND;
        while (Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_PEND) || Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_RUN)) {
            Thread.sleep(30 * 1000L);

            int numJobs = LibBat.lsb_openjobinfo(jobId, null, null, null, null, LibBat.ALL_JOB);
            try {
                Assert.assertEquals(numJobs, 1);
    
                IntByReference more = new IntByReference();

                jobInfoEnt jobInfo = LibBat.lsb_readjobinfo(more);
                Assert.assertNotNull(jobInfo, "Job info is null");
                Assert.assertEquals(more.getValue(), 0, "More job info results than expected");

                jobStatus = jobInfo.status;
            } finally {
                LibBat.lsb_closejobinfo();
            }
        }
        Assert.assertTrue(Utils.isFlagSet(jobStatus, LibBat.JOB_STAT_DONE), String.format("Unexpected job status: 0x%02x", jobStatus));

        Assert.assertTrue(FileUtils.waitFor(outFile, 120), "File not found: " + outFile.getAbsolutePath());
        System.out.println("--- output ---");
        System.out.println(FileUtils.readFileToString(outFile));
        System.out.println("--- output ---");
        Assert.assertTrue(outFile.delete(), "Unable to delete " + outFile.getAbsolutePath());
        Assert.assertEquals(reply.queue, req.queue, "LSF reply queue does not match requested queue.");
        System.out.println("Validating that we reached the end of the test without exit.");
    }
}
