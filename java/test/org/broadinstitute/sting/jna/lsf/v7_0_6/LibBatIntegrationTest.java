/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.jna.lsf.v7_0_6;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.jna.lsf.v7_0_6.LibBat.*;

import java.io.File;

/**
 * Really a unit test, but this test will only run on systems with LSF setup.
 */
public class LibBatIntegrationTest extends BaseTest {
    @Test
    public void testClusterName() {
        String clusterName = LibLsf.ls_getclustername();
        System.out.println("Cluster name: " + clusterName);
        Assert.assertNotNull(clusterName);
    }

    @Test
    public void testSubmitEcho() {
        String queue = "hour";
        File outFile = new File("LibBatIntegrationTest.out");

        Assert.assertFalse(LibBat.lsb_init("LibBatIntegrationTest") < 0, LibBat.lsb_sperror("lsb_init() failed"));

        submit req = new submit();

        for (int i = 0; i < LibLsf.LSF_RLIM_NLIMITS; i++)
            req.rLimits[i] = LibLsf.DEFAULT_RLIMIT;

        req.projectName = "LibBatIntegrationTest";
        req.options |= LibBat.SUB_PROJECT_NAME;

        req.queue = queue;
        req.options |= LibBat.SUB_QUEUE;

        req.outFile = outFile.getPath();
        req.options |= LibBat.SUB_OUT_FILE;

        req.command = "echo \"Hello world.\"";
        req.options2 |= LibBat.SUB2_BSUB_BLOCK;

        submitReply reply = new submitReply();
        long jobId = LibBat.lsb_submit(req, reply);

        Assert.assertFalse(jobId < 0, LibBat.lsb_sperror("Error dispatching"));
        Assert.assertTrue(FileUtils.waitFor(outFile, 120), "File not found: " + outFile.getAbsolutePath());
        Assert.assertTrue(outFile.delete(), "Unable to delete " + outFile.getAbsolutePath());
        Assert.assertEquals(reply.queue, req.queue, "LSF reply queue does not match requested queue.");
    }
}
