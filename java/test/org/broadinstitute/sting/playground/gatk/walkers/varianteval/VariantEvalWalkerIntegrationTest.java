package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.*;

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
        HashMap<String, String> md5 = new HashMap<String, String>();
        md5.put("", "d6b8c2d6c37d42d1ca2288799a8bd8e4");
        md5.put("-A", "0294b2e3915e88dfe2547e9db64ed1b3");

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
        for ( Map.Entry<String, String> e : md5.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                            " --rodBind eval,Variants," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                            " -T VariantEval" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -L 1:10,000,000-11,000,000" +
                            " --outerr %s" +
                            " --supressDateInformation " + e.getKey(),
                    1, // just one output file
                    Arrays.asList(e.getValue()));
            List<File> result = executeTest("testEvalVariantROD", spec).getFirst();
        }
    }

    @Test
    public void testEvalVariantRODConfSix() {
        List<String> md5 = new ArrayList<String>();
        md5.add("85cfefcac2dfb06545792605a3043a52");

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
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " --rodBind eval,Variants," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantEval" +
                        " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                        " -L 1:10,000,000-11,000,000" +
                        " --outerr %s" +
                        " --supressDateInformation" +
                        " -minPhredConfidenceScore 60",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalVariantRODConfSixty", spec).getFirst();
    }

     @Test
    public void testEvalVariantRODOutputViolations() {
        List<String> md5 = new ArrayList<String>();
        md5.add("e24732ffd95a78385a2c6986d1d3a359");

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
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " --rodBind eval,Variants," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
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
        md5.add("010d1c7ce773b39f3de1355eb9682e4d");
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
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " --rodBind eval,Variants," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
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
        md5.add("7d5a98c01051f96a684a383786da3d76");
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
                "-T VariantEval -R " + oneKGLocation + "reference/human_b36_both.fasta " +
                "-B eval,Variants," + validationDataLocation + "UMichVsBroad.venn.set1Only.calls " +
                "-D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -hc /humgen/gsa-scr1/GATK_Data/1KG_gffs/NA12878.1kg.gff " +
                "-G " +
                "--supressDateInformation " +
                "-L 1:1-10,000,000 " +
                "--outerr %s",
                1, // just one output file
                md5);
        List<File> result = executeTest("testEvalMarksGenotypingExample", spec).getFirst();
    }

    @Test
	public void testEvalRuntimeWithLotsOfIntervals() {
        List<String> md5 = new ArrayList<String>();
        md5.add("d11ea079fc1835514d392056a2c2a28d");
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T VariantEval -R " + oneKGLocation + "reference/human_b36_both.fasta " +
                "-B eval,Variants," + validationDataLocation + "NA12878.pilot_3.all.geli.calls " +
                "-D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod " +
                "--supressDateInformation " +
                "-L /humgen/gsa-scr1/GATK_Data/thousand_genomes_alpha_redesign.targets.b36.interval_list " +
                "--outerr %s",
                1, // just one output file                                                                                                                                                          
                md5);
        List<File> result = executeTest("testEvalRuntimeWithLotsOfIntervals", spec).getFirst();
    }

    @Test
    public void testVCFVariantEvals() {
        HashMap<String, String> md5 = new HashMap<String, String>();
        md5.put("", "3dda57ac7a9c8f3800726c9affb9d9bd");
        md5.put("-A", "d985e61fd0d7fc34c9c1a553e2881c67");
        md5.put("-A --includeFilteredRecords", "434c60986aa54c5fd07c22df1910ec44");
        md5.put("-A --sampleName NA12878", "aff844b88f71824a6cd3cce553325b17");
        md5.put("-A -vcfInfoSelector AF=0.50", "9ab9fa5d89cd6e3278d0d2b13cabbd51");

        for ( Map.Entry<String, String> e : md5.entrySet() ) {
            WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                    "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                            " --rodBind eval,VCF," + validationDataLocation + "NA12878.example1.vcf" +
                            " -T VariantEval" +
                            " --DBSNP /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod" +
                            " -hc /humgen/gsa-scr1/GATK_Data/1KG_gffs/NA12878.1kg.gff" +
                            " -G" +
                            " -L 1:1-10,000" +
                            " --outerr %s" +
                            " --supressDateInformation " + e.getKey(),
                    1, // just one output file
                    Arrays.asList(e.getValue()));
            List<File> result = executeTest("testVCFVariantEvals", spec).getFirst();
        }
    }


}

