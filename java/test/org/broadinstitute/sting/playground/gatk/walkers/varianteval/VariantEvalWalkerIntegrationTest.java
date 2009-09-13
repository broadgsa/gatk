package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author aaron
 *         <p/>
 *         Class VariantEvalWalkerTest
 *         <p/>
 *         test out the variant eval walker under different runtime conditions.
 */
public class VariantEvalWalkerIntegrationTest extends WalkerTest {

    @Test
    public void testEvalVariantROD() {
        List<String> md5 = new ArrayList<String>();
        md5.add("013d13d405f33b0fe4038f5aea087e7f");

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
                        " --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantEval" +
                        " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                        " -L 1:10,000,000-11,000,000" +
                        " --outerr %s" +
                        " --supressDateInformation",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalVariantROD", spec).getFirst();
    }

    @Test
    public void testEvalVariantRODConfSix() {
        List<String> md5 = new ArrayList<String>();
        md5.add("422d94ab60a3db42b0a5dc9c254b78f4");

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
         * --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls \
         * -minConfidenceScore 6
         */
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R /broad/1KG/reference/human_b36_both.fasta" +
                        " --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantEval" +
                        " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                        " -L 1:10,000,000-11,000,000" +
                        " --outerr %s" +
                        " --supressDateInformation" +
                        " -minConfidenceScore 6",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalVariantRODConfSix", spec).getFirst();
    }

     @Test
    public void testEvalVariantRODOutputViolations() {
        List<String> md5 = new ArrayList<String>();
        md5.add("8e613d57fd770ff04eb6531f583f4448");

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
         * --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls \
         * --includeViolations
         */
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R /broad/1KG/reference/human_b36_both.fasta" +
                        " --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantEval" +
                        " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                        " -L 1:10,000,000-11,000,000" +
                        " --outerr %s" +
                        " --supressDateInformation" +
                        " --includeViolations",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalVariantRODOutputViolations", spec).getFirst();
    }

    @Test
    public void testEvalGenotypeROD() {
        List<String> md5 = new ArrayList<String>();
        md5.add("439c2983a6e1250cd845c686a0e5a085");
        /**
         * the above MD5 was calculated after running the following command:
         *
         * java -jar ./dist/GenomeAnalysisTK.jar \
         * -R /broad/1KG/reference/human_b36_both.fasta \
         * -T VariantEval \
         * --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod \
         * -L 1:10,000,000-11,000,000 \
         * --outerr myVariantEval \
         * --supressDateInformation \
         * --evalContainsGenotypes \
         * --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls \
         * --rodBind hapmap-chip,GFF,/humgen/gsa-scr1/GATK_Data/1KG_gffs/NA12878.1kg.gff
         */

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R /broad/1KG/reference/human_b36_both.fasta" +
                        " --rodBind eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
                        " -T VariantEval" +
                        " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                        " -L 1:10,000,000-11,000,000" +
                        " --outerr %s" +
                        " --supressDateInformation" +
                        " --evalContainsGenotypes" +
                        " --rodBind hapmap-chip,GFF,/humgen/gsa-scr1/GATK_Data/1KG_gffs/NA12878.1kg.gff",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalGenotypeROD", spec).getFirst();
    }
}

