package org.broadinstitute.sting.utils.genotype.glf;

import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;

import org.testng.annotations.Test;


/**
 * 
 * @author aaron 
 * 
 * Class GLFRecordUnitTest
 *
 * Test out the basics of a GLFRecord
 */
public class GLFRecordUnitTest extends BaseTest {

    @Test
    public void testConstructingGLFRecord() {
        double likelihoods[] = new double[10];
        for (int i = 0; i < 10; i++) {
            likelihoods[i] = 10.0;
        }
        GLFRecord rec = new GLFSingleCall("1",'A',1,100,(short)200,likelihoods);

        Assert.assertTrue("1".equals(rec.contig));
        Assert.assertEquals(rec.getRefBase().toChar(), 'A');
        Assert.assertEquals(rec.getPosition(), 1);
        Assert.assertEquals(rec.getMinimumLikelihood(), 10);
        Assert.assertEquals(rec.getRmsMapQ(), 200);
        Assert.assertEquals(rec.getReadDepth(), 100);

    }

}
