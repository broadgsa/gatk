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

package org.broadinstitute.gatk.engine.phonehome;

import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.activeregion.ActivityProfileState;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.jets3t.service.S3Service;
import org.jets3t.service.S3ServiceException;
import org.jets3t.service.ServiceException;
import org.jets3t.service.model.S3Object;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

public class GATKRunReportUnitTest extends BaseTest {
    private final static boolean DEBUG = false;
    private static final long S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING = 30 * 1000;
    private static final String AWS_DOWNLOADER_CREDENTIALS_PROPERTIES_FILE = privateTestDir + "phonehome/awsDownloaderCredentials.properties";

    private Walker walker;
    private Exception exception;
    private GenomeAnalysisEngine engine;
    private String downloaderAccessKey;
    private String downloaderSecretKey;

    @BeforeClass
    public void setup() throws Exception {
        walker = new RunReportDummyReadWalker();
        exception = new IllegalArgumentException("javaException");
        engine = new GenomeAnalysisEngine();
        engine.setArguments(new GATKArgumentCollection());

        Properties awsProperties = new Properties();
        awsProperties.load(new FileInputStream(AWS_DOWNLOADER_CREDENTIALS_PROPERTIES_FILE));
        downloaderAccessKey = awsProperties.getProperty("accessKey");
        downloaderSecretKey = awsProperties.getProperty("secretKey");
    }

    @Test(enabled = ! DEBUG)
    public void testAWSKeysAreValid() {
        // throws an exception if they aren't
        GATKRunReport.checkAWSAreValid();
    }

    @Test(enabled = ! DEBUG)
    public void testAccessKey() throws Exception {
        testAWSKey(GATKRunReport.getAWSUploadAccessKey(), GATKRunReport.AWS_ACCESS_KEY_MD5);
    }

    @Test(enabled = ! DEBUG)
    public void testSecretKey() throws Exception {
        testAWSKey(GATKRunReport.getAWSUploadSecretKey(), GATKRunReport.AWS_SECRET_KEY_MD5);
    }

    private void testAWSKey(final String accessKey, final String expectedMD5) throws Exception {
        Assert.assertNotNull(accessKey, "AccessKey should not be null");
        final String actualmd5 = Utils.calcMD5(accessKey);
        Assert.assertEquals(actualmd5, expectedMD5);
    }

    @DataProvider(name = "GATKReportCreationTest")
    public Object[][] makeGATKReportCreationTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Walker readWalker = new RunReportDummyReadWalker();
        final Walker lociWalker = new RunReportDummyLocusWalker();
        final Walker rodWalker = new RunReportDummyRodWalker();
        final Walker artWalker = new RunReportDummyActiveRegionWalker();

        final Exception noException = null;
        final Exception javaException = new IllegalArgumentException("javaException");
        final Exception stingException = new ReviewedGATKException("GATKException");
        final Exception userException = new UserException("userException");

        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setArguments(new GATKArgumentCollection());

        for ( final Walker walker : Arrays.asList(readWalker, lociWalker, rodWalker, artWalker) ) {
            for ( final Exception exception : Arrays.asList(noException,  javaException, stingException, userException) ) {
                tests.add(new Object[]{walker, exception, engine});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "GATKReportCreationTest")
    public void testGATKReportCreationReadingAndWriting(final Walker walker, final Exception exception, final GenomeAnalysisEngine engine) throws Exception {
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.STDOUT);
        final ByteArrayOutputStream captureStream = new ByteArrayOutputStream();
        final boolean succeeded = report.postReportToStream(captureStream);
        Assert.assertTrue(succeeded, "Failed to write report to stream");
        Assert.assertFalse(report.exceptionOccurredDuringPost(), "Post succeeded but report says it failed");
        Assert.assertNull(report.getErrorMessage(), "Post succeeded but there was an error message");
        Assert.assertNull(report.getErrorThrown(), "Post succeeded but there was an error message");
        final InputStream readStream = new ByteArrayInputStream(captureStream.toByteArray());

        GATKRunReport deserialized = null;
        try {
            deserialized = GATKRunReport.deserializeReport(readStream);
        } catch ( Exception e ) {
            final String reportString = new String(captureStream.toByteArray());
            Assert.fail("Failed to deserialize GATK report " + reportString + " with exception " + e);
        }

        if ( deserialized != null )
            Assert.assertEquals(report, deserialized);
    }

    @DataProvider(name = "GATKAWSReportMode")
    public Object[][] makeGATKAWSReportMode() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final GATKRunReport.AWSMode mode : GATKRunReport.AWSMode.values() ) {
            tests.add(new Object[]{mode});
        }

        return tests.toArray(new Object[][]{});
    }

    // Will fail with timeout if AWS time out isn't working
    // Will fail with exception if AWS doesn't protect itself from errors
    @Test(enabled = ! DEBUG, dataProvider = "GATKAWSReportMode", timeOut = S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING * 2)
    public void testAWS(final GATKRunReport.AWSMode awsMode) {
        logger.warn("Starting testAWS mode=" + awsMode);

        // Use a shorter timeout than usual when we're testing GATKRunReport.AWSMode.TIMEOUT
        final long thisTestS3Timeout = awsMode == GATKRunReport.AWSMode.TIMEOUT ? 30 * 1000 : S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING;
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.AWS, thisTestS3Timeout);
        report.sendAWSToTestBucket();
        report.setAwsMode(awsMode);
        final S3Object s3Object = report.postReportToAWSS3();

        if ( awsMode == GATKRunReport.AWSMode.NORMAL ) {
            Assert.assertNotNull(s3Object, "Upload to AWS failed, s3Object was null. error was " + report.formatError());
            Assert.assertFalse(report.exceptionOccurredDuringPost(), "The upload should have succeeded but the report says it didn't.  Error was " + report.formatError());
            Assert.assertNull(report.getErrorMessage(), "Report succeeded but an error message was found");
            Assert.assertNull(report.getErrorThrown(), "Report succeeded but an thrown error was found");
            try {
                final GATKRunReport deserialized = GATKRunReport.deserializeReport(downloaderAccessKey, downloaderSecretKey, report.getS3ReportBucket(), s3Object);
                Assert.assertEquals(report, deserialized);
                deleteFromS3(report);
            } catch ( Exception e ) {
                Assert.fail("Failed to read, deserialize, or delete GATK report " + s3Object.getName() + " with exception " + e);
            }
        } else {
            Assert.assertNull(s3Object, "AWS upload should have failed for mode " + awsMode + " but got non-null s3 object back " + s3Object + " error was " + report.formatError());
            Assert.assertTrue(report.exceptionOccurredDuringPost(), "S3 object was null but the report says that the upload succeeded");
            Assert.assertNotNull(report.getErrorMessage(), "Report succeeded but an error message wasn't found");
            if ( awsMode == GATKRunReport.AWSMode.FAIL_WITH_EXCEPTION )
                Assert.assertNotNull(report.getErrorThrown());
        }
    }

    private void deleteFromS3(final GATKRunReport report) throws Exception {
        final S3Service s3Service = GATKRunReport.initializeAWSService(downloaderAccessKey, downloaderSecretKey);
        // Retrieve the whole data object we created previously
        s3Service.deleteObject(report.getS3ReportBucket(), report.getReportFileName());
    }

    @DataProvider(name = "PostReportByType")
    public Object[][] makePostReportByType() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final GATKRunReport.PhoneHomeOption et : GATKRunReport.PhoneHomeOption.values() ) {
            tests.add(new Object[]{et});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = ! DEBUG, dataProvider = "PostReportByType", timeOut = S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING * 2)
    public void testPostReportByType(final GATKRunReport.PhoneHomeOption type) {
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.AWS, S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING);
        Assert.assertFalse(report.exceptionOccurredDuringPost(), "An exception occurred during posting the report");
        final boolean succeeded = report.postReport(type);

        if ( type == GATKRunReport.PhoneHomeOption.NO_ET )
            Assert.assertFalse(succeeded, "NO_ET option shouldn't write a report");
        else {
            Assert.assertTrue(succeeded, "Any non NO_ET option should succeed in writing a report");

            if ( type == GATKRunReport.PhoneHomeOption.STDOUT ) {
                // nothing to do
            } else {
                // must have gone to AWS
                try {
                    Assert.assertTrue(report.wentToAWS(), "The report should have gone to AWS but the report says it wasn't");
                    deleteFromS3(report);
                } catch ( Exception e ) {
                    Assert.fail("Failed delete GATK report " + report.getReportFileName() + " with exception " + e);
                }
            }
        }
    }

    public interface S3Op {
        public void apply() throws ServiceException;
    }

    // Will fail with timeout if AWS time out isn't working
    // Will fail with exception if AWS doesn't protect itself from errors
    @Test(timeOut = S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING * 2)
    public void testAWSPublicKeyHasAccessControls() throws Exception {
        final GATKRunReport report = new GATKRunReport(walker, exception, engine, GATKRunReport.PhoneHomeOption.AWS, S3_PUT_TIMEOUT_IN_MILLISECONDS_FOR_TESTING);
        report.sendAWSToTestBucket();
        final S3Object s3Object = report.postReportToAWSS3();
        Assert.assertNotNull(s3Object, "Upload to AWS failed, s3Object was null. error was " + report.formatError());

        // create a service with the public key, and make sure it cannot list or delete
        final S3Service s3Service = GATKRunReport.initializeAWSService(GATKRunReport.getAWSUploadAccessKey(), GATKRunReport.getAWSUploadSecretKey());
        assertOperationNotAllowed("listAllBuckets", new S3Op() {
            @Override
            public void apply() throws S3ServiceException {
                s3Service.listAllBuckets();
            }
        });
        assertOperationNotAllowed("listBucket", new S3Op() {
            @Override
            public void apply() throws S3ServiceException { s3Service.listObjects(report.getS3ReportBucket()); }
        });
        assertOperationNotAllowed("createBucket", new S3Op() {
            @Override
            public void apply() throws S3ServiceException { s3Service.createBucket("ShouldNotCreate"); }
        });
        assertOperationNotAllowed("deleteObject", new S3Op() {
            @Override
            public void apply() throws ServiceException { s3Service.deleteObject(report.getS3ReportBucket(), report.getReportFileName()); }
        });
    }

    private void assertOperationNotAllowed(final String name, final S3Op op) {
        try {
            op.apply();
            // only gets here if the operation was successful
            Assert.fail("Operation " + name + " ran successfully but we expected to it fail");
        } catch ( ServiceException e ) {
            Assert.assertEquals(e.getErrorCode(), "AccessDenied");
        }
    }

    class RunReportDummyReadWalker extends ReadWalker<Integer, Integer> {
        @Override
        public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
            return 0;
        }

        @Override
        public Integer reduceInit() {
            return 0;
        }

        @Override
        public Integer reduce(Integer value, Integer sum) {
            return 0;
        }
    }

    class RunReportDummyLocusWalker extends LocusWalker<Integer, Integer> {
        @Override
        public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            return 0;
        }

        @Override
        public Integer reduceInit() {
            return 0;
        }

        @Override
        public Integer reduce(Integer value, Integer sum) {
            return 0;
        }
    }

    class RunReportDummyRodWalker extends RodWalker<Integer, Integer> {
        @Override
        public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            return 0;
        }

        @Override
        public Integer reduceInit() {
            return 0;
        }

        @Override
        public Integer reduce(Integer value, Integer sum) {
            return 0;
        }
    }

    class RunReportDummyActiveRegionWalker extends ActiveRegionWalker<Integer, Integer> {
        @Override
        public ActivityProfileState isActive(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            return new ActivityProfileState(ref.getLocus(), 0.0);
        }

        @Override
        public Integer map(ActiveRegion activeRegion, RefMetaDataTracker metaDataTracker) {
            return 0;
        }

        @Override
        public Integer reduceInit() {
            return 0;
        }

        @Override
        public Integer reduce(Integer value, Integer sum) {
            return 0;
        }
    }
}
