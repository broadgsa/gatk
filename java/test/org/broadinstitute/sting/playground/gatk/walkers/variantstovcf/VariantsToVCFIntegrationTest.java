package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.List;
import java.util.ArrayList;
import java.io.File;


/**
 * @author aaron
 *         <p/>
 *         Class VariantsToVCFIntegrationTest
 *         <p/>
 *         test(s) for the VariantsToVCF walker.
 */
public class VariantsToVCFIntegrationTest extends WalkerTest {


    @Test
    public void testVariantsToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("0b96a8046d2a06bd87f57df8bac1678d");

        /**
         * the above MD5 was calculated from running the following command:
         *
         * java -jar ./dist/GenomeAnalysisTK.jar \
         * -R /broad/1KG/reference/human_b36_both.fasta \
         * -T VariantEval \
         * --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod \
         * -L 1:10,000,000-11,000,000 \
         * --outerr myVariantEval \
         * --supressDateInformation \
         * --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls
         *
         */
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R /broad/1KG/reference/human_b36_both.fasta" +
                        " --rodBind NA123AB,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " --vcfout %s",
                1, // just one output file
                md5);
        List<File> result = executeTest("testVariantsToVCFUsingGeliInput", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("09660faa7cfad8af36602f79461c0605");

        /**
         * the above MD5 was calculated from running the following command:
         *
         * java -jar ./dist/GenomeAnalysisTK.jar \
         * -R /broad/1KG/reference/human_b36_both.fasta \
         * -T VariantEval \
         * --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod \
         * -L 1:10,000,000-11,000,000 \
         * --outerr myVariantEval \
         * --supressDateInformation \
         * --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls
         *
         */
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R /broad/1KG/reference/human_b36_both.fasta" +
                        " --rodBind NA123AB,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " --vcfout %s",
                1, // just one output file
                md5);
        List<File> result = executeTest("testVariantsToVCFUsingGeliInput", spec).getFirst();
    }

}
