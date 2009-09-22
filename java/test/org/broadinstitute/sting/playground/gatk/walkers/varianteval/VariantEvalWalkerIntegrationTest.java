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
        md5.add("7d1f8ddbf2d2f721a0c6dc39c97b584c");

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
        md5.add("27969a841874fc00d007b01ee008eb58");

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
        md5.add("63940f8d977b15bb1bcc9223ebeacf43");

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
        md5.add("c616348c1041157acd2715285c5b38b8");
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

    @Test
    public void testEvalMarksGenotypingExample() {
        List<String> md5 = new ArrayList<String>();
        md5.add("004cdf0c3ca46744daf42645785acc28");
        /**
         * Run with the following commands:
         * 
         * java -Xmx2048m -jar /humgen/gsa-hphome1/depristo/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar
         * -T VariantEval -R /broad/1KG/reference/human_b36_both.fasta -l INFO
         * -B eval,Variants,/humgen/gsa-scr1/ebanks/concordanceForMark/UMichVsBroad.venn.set1Only.calls
         * -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -hc /humgen/gsa-scr1/GATK_Data/1KG_gffs/NA12878.1kg.gff
         * -G -L 1 -o /humgen/gsa-scr1/ebanks/concordanceForMark/UMichVsBroad.venn.set1Only.calls.eval
         */

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T VariantEval -R /broad/1KG/reference/human_b36_both.fasta " +
                "-B eval,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/UMichVsBroad.venn.set1Only.calls " +
                "-D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -hc /humgen/gsa-scr1/GATK_Data/1KG_gffs/NA12878.1kg.gff " +
                "-G " +
                "--supressDateInformation " +
                "-L 1:1-10,000,000 " +
                "--outerr %s",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalMarksGenotypingExample", spec).getFirst();
    }
}

