package org.broadinstitute.sting.utils.genotype.glf;

import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Test;


/**
 * 
 * @author aaron 
 * 
 * Class GLFRecordTest
 *
 * Test out the basics of a GLFRecord
 */
public class GLFRecordTest extends BaseTest {

    @Test
    public void testConstructingGLFRecord() {
        double likelihoods[] = new double[10];
        for (int i = 0; i < 10; i++) {
            likelihoods[i] = 10.0;
        }
        GLFRecord rec = new GLFSingleCall("1",'A',1,100,(short)200,likelihoods);

        Assert.assertTrue("1".equals(rec.contig));
        Assert.assertEquals('A',rec.getRefBase().toChar());
        Assert.assertEquals(1,rec.getPosition());
        Assert.assertEquals(10,rec.getMinimumLikelihood());
        Assert.assertEquals(200,rec.getRmsMapQ());
        Assert.assertEquals(100,rec.getReadDepth());

    }

}
