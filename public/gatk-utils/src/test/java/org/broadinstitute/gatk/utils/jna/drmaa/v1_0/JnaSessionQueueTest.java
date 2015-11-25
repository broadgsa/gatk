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

package org.broadinstitute.gatk.utils.jna.drmaa.v1_0;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.utils.BaseTest;
import org.ggf.drmaa.*;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class JnaSessionQueueTest extends BaseTest {
    private String implementation = null;
    private static final SessionFactory factory = new JnaSessionFactory();

    @Test
    public void testDrmaa() throws Exception {
        Session session = factory.getSession();
        Version version = session.getVersion();
        System.out.println(String.format("DRMAA version: %d.%d", version.getMajor(), version.getMinor()));
        System.out.println(String.format("DRMAA contact(s): %s", session.getContact()));
        System.out.println(String.format("DRM system(s): %s", session.getDrmSystem()));
        System.out.println(String.format("DRMAA implementation(s): %s", session.getDrmaaImplementation()));
        this.implementation = session.getDrmaaImplementation();
    }

    @Test(dependsOnMethods = { "testDrmaa" })
    public void testSubmitEcho() throws Exception {
        if ( ! queueTestRunModeIsSet ) {
            throw new SkipException("Skipping testSubmitEcho because we are in queue test dry run mode");
        }

        if (implementation.contains("LSF")) {
            System.err.println("    ***********************************************************");
            System.err.println("   *************************************************************");
            System.err.println("   ****                                                     ****");
            System.err.println("  ****  Skipping JnaSessionQueueTest.testSubmitEcho()        ****");
            System.err.println("  ****      Are you using the dotkit .combined_LSF_SGE?      ****");
            System.err.println("   ****                                                     ****");
            System.err.println("   *************************************************************");
            System.err.println("    ***********************************************************");
            throw new SkipException("Skipping testSubmitEcho because correct DRMAA implementation not found");
        }

        File outFile = tryCreateNetworkTempFile("JnaSessionQueueTest.out");
        Session session = factory.getSession();
        session.init(null);
        try {
            JobTemplate template = session.createJobTemplate();
            template.setRemoteCommand("sh");
            template.setOutputPath(":" + outFile.getAbsolutePath());
            template.setJoinFiles(true);
            template.setArgs(Arrays.asList("-c", "echo \"Hello world.\""));

            String jobId = session.runJob(template);
            System.out.println(String.format("Job id %s", jobId));
            session.deleteJobTemplate(template);

            System.out.println("Waiting for job to run: " + jobId);
            int remotePs = Session.QUEUED_ACTIVE;

            List<Integer> runningStatuses = Arrays.asList(Session.QUEUED_ACTIVE, Session.RUNNING);

            while (runningStatuses.contains(remotePs)) {
                Thread.sleep(30 * 1000L);
                remotePs = session.getJobProgramStatus(jobId);
            }

            Assert.assertEquals(remotePs, Session.DONE, "Job status is not DONE.");

            JobInfo jobInfo = session.wait(jobId, Session.TIMEOUT_NO_WAIT);

            Assert.assertTrue(jobInfo.hasExited(), String.format("Job did not exit cleanly: %s", jobId));
            Assert.assertEquals(jobInfo.getExitStatus(), 0, String.format("Exit status for jobId %s is non-zero", jobId));
            if (jobInfo.hasSignaled())
                Assert.fail(String.format("JobId %s exited with signal %s and core dump flag %s", jobId, jobInfo.getTerminatingSignal(), jobInfo.hasCoreDump()));
            Assert.assertFalse(jobInfo.wasAborted(), String.format("Job was aborted: %s", jobId));
        } finally {
            session.exit();
        }

        Assert.assertTrue(FileUtils.waitFor(outFile, 120), "File not found: " + outFile.getAbsolutePath());
        System.out.println("--- output ---");
        System.out.println(FileUtils.readFileToString(outFile));
        System.out.println("--- output ---");
        Assert.assertTrue(outFile.delete(), "Unable to delete " + outFile.getAbsolutePath());
        System.out.println("Validating that we reached the end of the test without exit.");
    }

    @Test
    public void testCollectionConversions() {
        Collection<String> list = Arrays.asList("a=1", "foo=bar", "empty=");
        Map<String, String> map = new LinkedHashMap<String, String>();
        map.put("a", "1");
        map.put("foo", "bar");
        map.put("empty", "");

        Assert.assertEquals(JnaSession.collectionToMap(list), map);
        Assert.assertEquals(JnaSession.mapToCollection(map), list);
    }

    @Test
    public void testLimitConversions() {
        Assert.assertEquals(JnaSession.formatLimit(0), "0:00:00");
        Assert.assertEquals(JnaSession.formatLimit(59), "0:00:59");
        Assert.assertEquals(JnaSession.formatLimit(60), "0:01:00");
        Assert.assertEquals(JnaSession.formatLimit(3540), "0:59:00");
        Assert.assertEquals(JnaSession.formatLimit(3599), "0:59:59");
        Assert.assertEquals(JnaSession.formatLimit(7200), "2:00:00");
        Assert.assertEquals(JnaSession.formatLimit(7260), "2:01:00");
        Assert.assertEquals(JnaSession.formatLimit(7261), "2:01:01");

        Assert.assertEquals(JnaSession.parseLimit("0"), 0);
        Assert.assertEquals(JnaSession.parseLimit("00"), 0);
        Assert.assertEquals(JnaSession.parseLimit("0:00"), 0);
        Assert.assertEquals(JnaSession.parseLimit("00:00"), 0);
        Assert.assertEquals(JnaSession.parseLimit("0:00:00"), 0);

        Assert.assertEquals(JnaSession.parseLimit("1"), 1);
        Assert.assertEquals(JnaSession.parseLimit("01"), 1);
        Assert.assertEquals(JnaSession.parseLimit("0:01"), 1);
        Assert.assertEquals(JnaSession.parseLimit("00:01"), 1);
        Assert.assertEquals(JnaSession.parseLimit("0:00:01"), 1);

        Assert.assertEquals(JnaSession.parseLimit("10"), 10);
        Assert.assertEquals(JnaSession.parseLimit("0:10"), 10);
        Assert.assertEquals(JnaSession.parseLimit("00:10"), 10);
        Assert.assertEquals(JnaSession.parseLimit("0:00:10"), 10);

        Assert.assertEquals(JnaSession.parseLimit("1:0"), 60);
        Assert.assertEquals(JnaSession.parseLimit("1:00"), 60);
        Assert.assertEquals(JnaSession.parseLimit("01:00"), 60);
        Assert.assertEquals(JnaSession.parseLimit("0:01:00"), 60);

        Assert.assertEquals(JnaSession.parseLimit("1:00:00"), 3600);

        Assert.assertEquals(JnaSession.parseLimit("1:02:03"), 3723);
    }
}
