package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Nov 17, 2009
 * Time: 6:58:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class SecondBaseSkewIntegrationTest extends WalkerTest {

    @Test
    public void secondBaseSkewTest() {
        for ( int test = 1; test <= 5; test ++ ) {
            String bamFilePath = VariantAnnotatorIntegrationTest.validationDataPath()+VariantAnnotatorIntegrationTest.secondBaseTestFile(test)+".a2b.recal.annotation_subset.bam";
            String callFile = VariantAnnotatorIntegrationTest.validationDataPath()+VariantAnnotatorIntegrationTest.secondBaseTestFile(test)+".a2b.ssg1b.geli.calls";
            String args = VariantAnnotatorIntegrationTest.secondBaseTestString()+" -I "+bamFilePath+" -B variant,Variants,"+callFile+" "+VariantAnnotatorIntegrationTest.secondBaseTestInterval(test)+" -sample variant";
            WalkerTestSpec spec = new WalkerTestSpec(args,1, Arrays.asList(VariantAnnotatorIntegrationTest.secondBaseTestmd5(test)));
            executeTest("Second base skew annotator, test number "+Integer.toString(test), spec);
        }
    }

}
