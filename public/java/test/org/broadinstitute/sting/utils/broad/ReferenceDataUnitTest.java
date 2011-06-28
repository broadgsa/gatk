package org.broadinstitute.sting.utils.broad;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class ReferenceDataUnitTest {
    @Test
    public void testNames() {
        Assert.assertEquals(ReferenceData.HG18.getName(), "hg18");
        Assert.assertEquals(ReferenceData.HG19.getName(), "hg19");
    }

    @Test
    public void testFilesExist() {
        for (ReferenceData data: ReferenceData.values()) {
            Assert.assertTrue(new File(data.getReference()).exists());
            Assert.assertTrue(new File(data.getRefseq()).exists());
            for (int version: data.getDbsnpVersions()) {
                Assert.assertTrue(new File(data.getDbsnp(version)).exists());
            }
        }
    }

    @Test
    public void testDbsnps() {
        Assert.assertTrue(new File(ReferenceData.HG18.getDbsnp(129)).exists());
        Assert.assertTrue(new File(ReferenceData.HG19.getDbsnp(129)).exists());
        Assert.assertTrue(new File(ReferenceData.HG19.getDbsnp(132)).exists());
        Assert.assertNull(ReferenceData.HG19.getDbsnp(130));
    }

    @Test
    public void testDbsnpTypes() {
        Assert.assertEquals(ReferenceData.HG18.getDbsnpType(129), "DBSNP");
        Assert.assertEquals(ReferenceData.HG19.getDbsnpType(129), "VCF");
        Assert.assertEquals(ReferenceData.HG19.getDbsnpType(132), "VCF");
        Assert.assertNull(ReferenceData.HG19.getDbsnpType(130));
    }

    @Test
    public void testGetByReference() {
        Assert.assertEquals(ReferenceData.getByReference(BaseTest.hg18Reference), ReferenceData.HG18);
        Assert.assertEquals(ReferenceData.getByReference(BaseTest.hg19Reference), ReferenceData.HG19);
        Assert.assertEquals(ReferenceData.getByReference("none"), null);
    }
}
