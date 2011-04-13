package org.broadinstitute.sting.utils.broad;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;

public class PicardAggregationUtilsUnitTest {
    public static final String PROJECT = "C474";
    public static final String SAMPLE = "NA19651";
    public static final String MISSING_PROJECT = "C0";
    public static final String MISSING_SAMPLE = "0";
    private int latestVersion = -1;

    @Test
    public void testGetLatestVersion() {
        latestVersion = PicardAggregationUtils.getLatestVersion(PROJECT, SAMPLE);
        System.out.println(String.format("Latest version for %s %s is %d", PROJECT, SAMPLE, latestVersion));
        Assert.assertTrue(latestVersion > 0);
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(PROJECT, SAMPLE, latestVersion), latestVersion);
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testGetSampleBam() throws Exception {
        String test = PicardAggregationUtils.getSampleBam(PROJECT, SAMPLE);
        String latest = PicardAggregationUtils.getSampleBam(PROJECT, SAMPLE, latestVersion);
        Assert.assertEquals(test, latest);
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testGetSampleDir() throws Exception {
        String test = PicardAggregationUtils.getSampleDir(PROJECT, SAMPLE);
        String latest = PicardAggregationUtils.getSampleDir(PROJECT, SAMPLE, latestVersion);
        Assert.assertEquals(test, latest);
    }

    @Test(dependsOnMethods = "testGetLatestVersion")
    public void testIsFinished() {
        Assert.assertTrue(PicardAggregationUtils.isFinished(PROJECT, SAMPLE, latestVersion));
        Assert.assertFalse(PicardAggregationUtils.isFinished(PROJECT, SAMPLE, latestVersion + 1));
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testMissingSampleBam() throws Exception {
        PicardAggregationUtils.getSampleBam(MISSING_PROJECT, MISSING_SAMPLE);
    }

    @Test(expectedExceptions = FileNotFoundException.class)
    public void testMissingSampleDir() throws Exception {
        PicardAggregationUtils.getSampleDir(MISSING_PROJECT, MISSING_SAMPLE);
    }

    @Test
    public void testLatestVersionMissing() {
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(MISSING_PROJECT, MISSING_SAMPLE), 0);
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(MISSING_PROJECT, MISSING_SAMPLE, -1), -1);
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(MISSING_PROJECT, MISSING_SAMPLE, 0), 0);
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(MISSING_PROJECT, MISSING_SAMPLE, 1), 1);
        Assert.assertEquals(PicardAggregationUtils.getLatestVersion(MISSING_PROJECT, MISSING_SAMPLE, 2), 2);
    }
}
