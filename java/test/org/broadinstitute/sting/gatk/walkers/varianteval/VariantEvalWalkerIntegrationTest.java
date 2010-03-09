package org.broadinstitute.sting.gatk.walkers.varianteval;

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
        md5.put("", "9cfda40f521d75a3e8bafc44a663c14a");
        md5.put("-A", "8fea7cc25f551ce170636fc35c5ae0fe");

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
        md5.add("11d636d105f902680c46b9f2e330d922");

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
        md5.add("12ecd457d62329e9d4e593de904a457d");

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
        md5.add("6ed44fd586c89dafd40cb8e0194dc456");
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
        md5.add("c0396cfe89a63948aebbbae0a0e06678");
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
        md5.add("6a90341517fc3c5026529301d9970c7b");
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
        md5.put("", "ee6b096169d6c5e2ce49d394fbec799b");
        md5.put("-A", "a443193c0810363f85278b1cfaed2fff");
        md5.put("-A --includeFilteredRecords", "812d7f2ecac28b1be7e7028af17df9c0");
        md5.put("-A --sampleName NA12878", "a443193c0810363f85278b1cfaed2fff");
        md5.put("-A -vcfInfoSelector AF=0.50", "afed4bf0c9f11b86f6e5356012f9cf2d");

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

