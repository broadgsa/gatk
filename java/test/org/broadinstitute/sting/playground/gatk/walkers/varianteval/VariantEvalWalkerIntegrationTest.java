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
        md5.put("", "4b5a43cb4fa1b82e1bb361632754eaa1");
        md5.put("-A", "7dd6aa5379ea5422fe7db4e79baad1f7");

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
        md5.add("cfcc2a8e0587cb288f0472118a9bdc1d");

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
        md5.add("ce708258676bbd4e71e8714828e7d695");

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
        md5.add("4df6a9c9d09a7bf034f29db24fb7a923");
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
        md5.add("f79448f380cfb56be3585547b794a48c");
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
        md5.add("4d01485ed2d3cb5edbff24cc0c05b205");
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
        md5.put("", "9702cba801cf216324607c8ecb7871e4");
        md5.put("-A", "e9e02eba2e3423bde18ddb0a33a50e34");
        md5.put("-A --includeFilteredRecords", "f0cb51e2c7917e502ca510cc9a3382a9");
        md5.put("-A --sampleName NA12878", "b936723b0de69ed2698ab2af4b70a1e8");
        md5.put("-A -vcfInfoSelector AF=0.50", "a8299c6ba370bf9334327ec03e8885b7");

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

