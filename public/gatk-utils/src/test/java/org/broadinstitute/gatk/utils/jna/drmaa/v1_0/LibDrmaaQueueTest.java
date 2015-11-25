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

import com.sun.jna.Memory;
import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;
import com.sun.jna.StringArray;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class LibDrmaaQueueTest extends BaseTest {
    private String implementation = null;

    @Test
    public void testDrmaa() throws Exception {
        Memory error = new Memory(LibDrmaa.DRMAA_ERROR_STRING_BUFFER);
        int errnum;

        IntByReference major = new IntByReference();
        IntByReference minor = new IntByReference();
        Memory contact = new Memory(LibDrmaa.DRMAA_CONTACT_BUFFER);
        Memory drmSystem = new Memory(LibDrmaa.DRMAA_DRM_SYSTEM_BUFFER);
        Memory drmaaImplementation = new Memory(LibDrmaa.DRMAA_DRMAA_IMPLEMENTATION_BUFFER);

        errnum = LibDrmaa.drmaa_version(major, minor, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

        if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
            Assert.fail(String.format("Could not get version from the DRMAA library: %s", error.getString(0)));

        System.out.println(String.format("DRMAA version: %d.%d", major.getValue(), minor.getValue()));

        errnum = LibDrmaa.drmaa_get_contact(contact, LibDrmaa.DRMAA_CONTACT_BUFFER_LEN, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

        if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
            Assert.fail(String.format("Could not get contacts from the DRMAA library: %s", error.getString(0)));

        System.out.println(String.format("DRMAA contact(s): %s", contact.getString(0)));

        errnum = LibDrmaa.drmaa_get_DRM_system(drmSystem, LibDrmaa.DRMAA_DRM_SYSTEM_BUFFER_LEN, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

        if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
            Assert.fail(String.format("Could not get DRM system from the DRMAA library: %s", error.getString(0)));

        System.out.println(String.format("DRM system(s): %s", drmSystem.getString(0)));

        errnum = LibDrmaa.drmaa_get_DRMAA_implementation(drmaaImplementation, LibDrmaa.DRMAA_DRMAA_IMPLEMENTATION_BUFFER_LEN, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

        if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
            Assert.fail(String.format("Could not get DRMAA implementation from the DRMAA library: %s", error.getString(0)));

        System.out.println(String.format("DRMAA implementation(s): %s", drmaaImplementation.getString(0)));

        this.implementation = drmaaImplementation.getString(0);
    }

    @Test(dependsOnMethods = { "testDrmaa" })
    public void testSubmitEcho() throws Exception {
        if ( ! queueTestRunModeIsSet ) {
            throw new SkipException("Skipping testSubmitEcho because we are in pipeline test dry run mode");
        }

        if (implementation.contains("LSF")) {
            System.err.println("    *********************************************************");
            System.err.println("   ***********************************************************");
            System.err.println("   ****                                                   ****");
            System.err.println("  ****  Skipping LibDrmaaQueueTest.testSubmitEcho()        ****");
            System.err.println("  ****     Are you using the dotkit .combined_LSF_SGE?     ****");
            System.err.println("   ****                                                   ****");
            System.err.println("   ***********************************************************");
            System.err.println("    *********************************************************");
            throw new SkipException("Skipping testSubmitEcho because correct DRMAA implementation not found");
        }

        Memory error = new Memory(LibDrmaa.DRMAA_ERROR_STRING_BUFFER);
        int errnum;

    File outFile = tryCreateNetworkTempFile("LibDrmaaQueueTest.out");

        errnum = LibDrmaa.drmaa_init(null, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

        if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
            Assert.fail(String.format("Could not initialize the DRMAA library: %s", error.getString(0)));

        try {
            PointerByReference jtRef = new PointerByReference();
            Pointer jt;
            Memory jobIdMem = new Memory(LibDrmaa.DRMAA_JOBNAME_BUFFER);
            String jobId;
            IntByReference remotePs = new IntByReference();
            IntByReference stat = new IntByReference();
            PointerByReference rusage = new PointerByReference();

            errnum = LibDrmaa.drmaa_allocate_job_template(jtRef, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Could not create job template: %s", error.getString(0)));

            jt = jtRef.getValue();

            errnum = LibDrmaa.drmaa_set_attribute(jt, LibDrmaa.DRMAA_REMOTE_COMMAND, "sh", error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Could not set attribute \"%s\": %s", LibDrmaa.DRMAA_REMOTE_COMMAND, error.getString(0)));

            errnum = LibDrmaa.drmaa_set_attribute(jt, LibDrmaa.DRMAA_OUTPUT_PATH, ":" + outFile.getAbsolutePath(), error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Could not set attribute \"%s\": %s", LibDrmaa.DRMAA_OUTPUT_PATH, error.getString(0)));

            errnum = LibDrmaa.drmaa_set_attribute(jt, LibDrmaa.DRMAA_JOIN_FILES, "y", error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Could not set attribute \"%s\": %s", LibDrmaa.DRMAA_JOIN_FILES, error.getString(0)));

            StringArray args = new StringArray(new String[] { "-c", "echo \"Hello world.\"" });

            errnum = LibDrmaa.drmaa_set_vector_attribute(jt, LibDrmaa.DRMAA_V_ARGV, args, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Could not set attribute \"%s\": %s", LibDrmaa.DRMAA_V_ARGV, error.getString(0)));

            errnum = LibDrmaa.drmaa_run_job(jobIdMem, LibDrmaa.DRMAA_JOBNAME_BUFFER_LEN, jt, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Could not submit job: %s", error.getString(0)));

            jobId = jobIdMem.getString(0);

            System.out.println(String.format("Job id %s", jobId));

            errnum = LibDrmaa.drmaa_delete_job_template(jt, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Could not delete job template: %s", error.getString(0)));

            System.out.println("Waiting for job to run: " + jobId);
            remotePs.setValue(LibDrmaa.DRMAA_PS.DRMAA_PS_QUEUED_ACTIVE);

            List<Integer> runningStatuses = Arrays.asList(
                    LibDrmaa.DRMAA_PS.DRMAA_PS_QUEUED_ACTIVE, LibDrmaa.DRMAA_PS.DRMAA_PS_RUNNING);

            while (runningStatuses.contains(remotePs.getValue())) {
                Thread.sleep(30 * 1000L);

                errnum = LibDrmaa.drmaa_job_ps(jobId, remotePs, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

                if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                    Assert.fail(String.format("Could not get status for jobId %s: %s", jobId, error.getString(0)));
            }

            Assert.assertEquals(remotePs.getValue(), LibDrmaa.DRMAA_PS.DRMAA_PS_DONE, "Job status is not DONE.");

            errnum = LibDrmaa.drmaa_wait(jobId, Pointer.NULL, new NativeLong(0), stat, LibDrmaa.DRMAA_TIMEOUT_NO_WAIT,
                    rusage, error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Wait failed for jobId %s: %s", jobId, error.getString(0)));

            IntByReference exited = new IntByReference();
            IntByReference exitStatus = new IntByReference();
            IntByReference signaled = new IntByReference();
            Memory signal = new Memory(LibDrmaa.DRMAA_SIGNAL_BUFFER);
            IntByReference coreDumped = new IntByReference();
            IntByReference aborted = new IntByReference();

            errnum = LibDrmaa.drmaa_wifexited(exited, stat.getValue(), error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Exit check failed for jobId %s: %s", jobId, error.getString(0)));

            Assert.assertTrue(exited.getValue() != 0, String.format("Job did not exit cleanly: %s", jobId));

            errnum = LibDrmaa.drmaa_wexitstatus(exitStatus, stat.getValue(), error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Exit status failed for jobId %s: %s", jobId, error.getString(0)));

            Assert.assertEquals(exitStatus.getValue(), 0, String.format("Exit status for jobId %s is non-zero", jobId));

            errnum = LibDrmaa.drmaa_wifsignaled(signaled, stat.getValue(), error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Signaled check failed for jobId %s: %s", jobId, error.getString(0)));

            if (signaled.getValue() != 0) {
                errnum = LibDrmaa.drmaa_wtermsig(signal, LibDrmaa.DRMAA_SIGNAL_BUFFER_LEN, stat.getValue(), error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

                if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                    Assert.fail(String.format("Signal lookup failed for jobId %s: %s", jobId, error.getString(0)));

                errnum = LibDrmaa.drmaa_wcoredump(coreDumped, stat.getValue(), error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

                if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                    Assert.fail(String.format("Core dump check failed for jobId %s: %s", jobId, error.getString(0)));

                Assert.fail(String.format("JobId %s exited with signal %s and core dump flag %d", jobId, signal.getString(0), coreDumped.getValue()));
            }

            errnum = LibDrmaa.drmaa_wifaborted(aborted, stat.getValue(), error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                Assert.fail(String.format("Aborted check failed for jobId %s: %s", jobId, error.getString(0)));

            Assert.assertTrue(aborted.getValue() == 0, String.format("Job was aborted: %s", jobId));

        } finally {
            if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS) {
                LibDrmaa.drmaa_exit(error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);
            } else {
                errnum = LibDrmaa.drmaa_exit(error, LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

                if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
                    Assert.fail(String.format("Could not shut down the DRMAA library: %s", error.getString(0)));
            }
        }

        Assert.assertTrue(FileUtils.waitFor(outFile, 120), "File not found: " + outFile.getAbsolutePath());
        System.out.println("--- output ---");
        System.out.println(FileUtils.readFileToString(outFile));
        System.out.println("--- output ---");
        Assert.assertTrue(outFile.delete(), "Unable to delete " + outFile.getAbsolutePath());
        System.out.println("Validating that we reached the end of the test without exit.");
    }
}
