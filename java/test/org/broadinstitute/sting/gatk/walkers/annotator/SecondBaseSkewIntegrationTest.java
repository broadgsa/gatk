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

    // TODO -- reinstitute these integration tests when GELI -> VariantContext is supported

    //@Test
    public void secondBaseSkewTest() {
        for ( int test = 1; test <= 2; test ++ ) {
            String bamFilePath = VariantAnnotatorIntegrationTest.validationDataPath()+VariantAnnotatorIntegrationTest.secondBaseTestFile(test)+".a2b.recal.annotation_subset.bam";
            String callFile = VariantAnnotatorIntegrationTest.validationDataPath()+VariantAnnotatorIntegrationTest.secondBaseTestFile(test)+".a2b.ssg1b.geli.calls";
            String args = VariantAnnotatorIntegrationTest.secondBaseTestString()+" -I "+bamFilePath+" -B variant,Variants,"+callFile+" "+VariantAnnotatorIntegrationTest.secondBaseTestInterval(test)+" -sample variant";
            WalkerTestSpec spec = new WalkerTestSpec(args,1, Arrays.asList(VariantAnnotatorIntegrationTest.secondBaseTestmd5(test)));
            executeTest("Second base skew annotator, test number "+Integer.toString(test), spec);
        }
    }

    //@Test
    public void testOnE2File() {
        String test_args = "-T VariantAnnotator -A SecondBaseSkew "
                +"-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta "
                +"-I " + validationDataLocation + "FHSP_pool3_2bannot.bam "
                +"-B variant,Variants," + validationDataLocation + "FHS_pilot_pool3_raw_calls.geli "
                +"-o %s -sample variant -L " + validationDataLocation + "FHS_test_intervals.interval_list";

        String md5_for_this_test = "c1116b3196cc9f553ae2442a4063bc5e";

        WalkerTestSpec spec = new WalkerTestSpec(test_args,1, Arrays.asList(md5_for_this_test));
        executeTest("Testing on E2 annotated but not Q2 annotated file ",spec);


    }

    //@Test
    public void testOnUnannotatedFile() {
        String test_args = "-T VariantAnnotator -A SecondBaseSkew "
                +"-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta "
                +"-I " + validationDataLocation + "FHSP_pool3_test.bam "
                +"-B variant,Variants," + validationDataLocation + "FHS_pilot_pool3_raw_calls.geli "
                +"-o %s -sample variant -L " + validationDataLocation + "FHS_test_intervals.interval_list";

        String md5_for_this_test = "a297259694ac88f769a45bce96a08e51";

        WalkerTestSpec spec = new WalkerTestSpec(test_args,1, Arrays.asList(md5_for_this_test));
        executeTest("Testing on bam file without 2bb annotations ",spec);
    }

    @Test
    public void testOnIndels() {
        String test_args = "-T VariantAnnotator -I " + validationDataLocation + "FHS_Pileup_Test.bam"
                     + " -R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -A SecondBaseSkew"
                     + " -sample variant -B variant,VCF," + validationDataLocation + "FHS_pileup_test_chr15.vcf"
                     + " -o %s -L chr15:46347148";
        String expected_md5 = "f4aff94a713fcf4dc5750d6a6a884ad3";
        WalkerTestSpec spec = new WalkerTestSpec(test_args,1,Arrays.asList(expected_md5));
        executeTest("Testing on locus with many indels", spec);
    }

    // todo -- chris needs to fix this
//    @Test
//    public void testPrimaryBaseSecondaryBaseOnIndels() {
//        String test_args = "-T VariantAnnotator -I /humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_Pileup_Test.bam"
//                     + " -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -A PrimaryBaseSecondaryBaseSymmetry"
//                     + " -sample variant -B variant,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/FHS_pileup_test_chr15.vcf"
//                     + " -vcf %s -L chr15:46347148";
//        String expected_md5 = "9b587be7a270c6df7e0affcfc61a861a";
//        WalkerTestSpec spec = new WalkerTestSpec(test_args,1,Arrays.asList(expected_md5));
//        executeTest("Testing PrimaryBaseSeHcondaryBaseSymmetry on locus with many indels", spec);
//    }

}
