package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {
        String[] md5DoC = {"c0a7e2fc07d565e633b3064f9f3cdaf5", "21c8e1f9dc65fdfb39347547f9b04011"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T VariantFiltration -X DepthOfCoverage:max=70 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5DoC));
        executeTest("testDoCFilter", spec1);

        String[] md5AlleleBalance = {"aa0f7800cfd346236620ae0eac220817", "a13e4ce6260bf9f33ca99dc808b8e6ad"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T VariantFiltration -X AlleleBalance:low=0.25,high=0.75 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5AlleleBalance));
        executeTest("testAlleleBalanceFilter", spec2);

        String[] md5Strand = {"9f430f251dbeb58a2f80a1306a5dd492", "0f7db0aad764268ee8fa3b857df8d87d"};
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T VariantFiltration -X FisherStrand:pvalue=0.0001 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Strand));
        executeTest("testStrandFilter", spec3);

        String[] md5Lod = {"56177258c0b3944c043f86faee4b42ae", "7e0c4f2b0fda85fd2891eee76c396a55"};
        WalkerTestSpec spec4 = new WalkerTestSpec(
                "-T VariantFiltration -X LodThreshold:lod=10 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Lod));
        executeTest("testLodFilter", spec4);

        String[] md5MQ0 = {"0e303c32f5c1503f4c875771f28fc46c", "3203de335621851bccf596242b079e23"};
        WalkerTestSpec spec5 = new WalkerTestSpec(
                "-T VariantFiltration -X MappingQualityZero:max=70 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5MQ0));
        executeTest("testMappingQuality0Filter", spec5);

        String[] md5MQ = {"946462a6199e9453784e0942e18e6830", "ecc777feedea61f7b570d114c2ab89b1"};
        WalkerTestSpec spec6 = new WalkerTestSpec(
                "-T VariantFiltration -X MappingQuality:min=20 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5MQ));
        executeTest("testRMSMappingQualityFilter", spec6);

        String[] md5OnOff = {"2ff84e104ce73e347e55d272170b4d03", "67f2e1bc025833b0fa31f47195198997"};
        WalkerTestSpec spec7 = new WalkerTestSpec(
                "-T VariantFiltration -X OnOffGenotypeRatio:threshold=0.9 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5OnOff));
        executeTest("testOnOffGenotypeFilter", spec7);

        String[] md5Clusters = {"e6a1c088678b1c31ff340ebd622b476e", "8fa6b6ffc93ee7fb8d6b52a7fb7815ef"};
        WalkerTestSpec spec8 = new WalkerTestSpec(
                "-T VariantFiltration -X ClusteredSnps:window=10,snps=3 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Clusters));
        executeTest("testClusteredSnpsFilter", spec8);

        String[] md5Indels = {"82e555b76c12474154f8e5e402516d73", "8e0e915a1cb63d7049e0671ed00101fe"};
        WalkerTestSpec spec9 = new WalkerTestSpec(
                "-T VariantFiltration -X IndelArtifact -B indels,PointIndel,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.indels -B cleaned,CleanedOutSNP,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.realigner_badsnps -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Indels));
        executeTest("testIndelArtifactFilter", spec9);






    }
}