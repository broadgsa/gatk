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
        for ( int test = 1; test <= 2; test ++ ) {
            String bamFilePath = VariantAnnotatorIntegrationTest.validationDataPath()+VariantAnnotatorIntegrationTest.secondBaseTestFile(test)+".a2b.recal.annotation_subset.bam";
            String callFile = VariantAnnotatorIntegrationTest.validationDataPath()+VariantAnnotatorIntegrationTest.secondBaseTestFile(test)+".a2b.ssg1b.geli.calls";
            String args = VariantAnnotatorIntegrationTest.secondBaseTestString()+" -I "+bamFilePath+" -B variant,Variants,"+callFile+" "+VariantAnnotatorIntegrationTest.secondBaseTestInterval(test)+" -sample variant";
            WalkerTestSpec spec = new WalkerTestSpec(args,1, Arrays.asList(VariantAnnotatorIntegrationTest.secondBaseTestmd5(test)));
            executeTest("Second base skew annotator, test number "+Integer.toString(test), spec);
        }
    }

    @Test
    public void testOnE2File() {
        String test_args = "-T VariantAnnotator -A SecondBaseSkew "
                +"-R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta "
                +"-I /humgen/gsa-scr1/GATK_Data/Validation_Data/FHSP_pool3_2bannot.bam "
                +"-B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_raw_calls.geli "
                +"-vcf %s -sample variant -L /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_test_intervals.interval_list";

        String md5_for_this_test = "6a71095e1cfade39d909b35f6c99d1ca";

        WalkerTestSpec spec = new WalkerTestSpec(test_args,1, Arrays.asList(md5_for_this_test));
        executeTest("Testing on E2 annotated but not Q2 annotated file ",spec);


    }

    @Test
    public void testOnUnannotatedFile() {
        String test_args = "-T VariantAnnotator -A SecondBaseSkew "
                +"-R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta "
                +"-I /humgen/gsa-scr1/GATK_Data/Validation_Data/FHSP_pool3_test.bam "
                +"-B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pilot_pool3_raw_calls.geli "
                +"-vcf %s -sample variant -L /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_test_intervals.interval_list";

        String md5_for_this_test = "f105fd8a7ae7026a55107b86e768553a";

        WalkerTestSpec spec = new WalkerTestSpec(test_args,1, Arrays.asList(md5_for_this_test));
        executeTest("Testing on bam file without 2bb annotations ",spec);
    }

}
