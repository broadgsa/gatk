package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        String[] md5DoC = {"5f176d27eba04c8e08b54b60fd96758e", "21c8e1f9dc65fdfb39347547f9b04011"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T VariantFiltration -X DepthOfCoverage:max=70 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5DoC));
        executeTest("testDoCFilter", spec1);

        String[] md5AlleleBalance = {"dfa47f9465214445b737ec49025d6b5f", "a13e4ce6260bf9f33ca99dc808b8e6ad"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T VariantFiltration -X AlleleBalance:low=0.25,high=0.75 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5AlleleBalance));
        executeTest("testAlleleBalanceFilter", spec2);

        String[] md5Strand = {"2d32ddd6daccd08f3df148f4d373fa5e", "0f7db0aad764268ee8fa3b857df8d87d"};
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T VariantFiltration -X FisherStrand:pvalue=0.0001 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Strand));
        executeTest("testStrandFilter", spec3);

        String[] md5Lod = {"157e7beab140f4252ed41ebfc9a508f2", "7e0c4f2b0fda85fd2891eee76c396a55"};
        WalkerTestSpec spec4 = new WalkerTestSpec(
                "-T VariantFiltration -X LodThreshold:lod=10 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Lod));
        executeTest("testLodFilter", spec4);

        String[] md5MQ0 = {"15738f3717dfa3192566a7ddf301fe11", "3203de335621851bccf596242b079e23"};
        WalkerTestSpec spec5 = new WalkerTestSpec(
                "-T VariantFiltration -X MappingQualityZero:max=70 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5MQ0));
        executeTest("testMappingQuality0Filter", spec5);

        String[] md5MQ = {"e7fb1f166d2359f653f98f8215215aca", "ecc777feedea61f7b570d114c2ab89b1"};
        WalkerTestSpec spec6 = new WalkerTestSpec(
                "-T VariantFiltration -X MappingQuality:min=20 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5MQ));
        executeTest("testRMSMappingQualityFilter", spec6);

        String[] md5OnOff = {"fe4acc1505f19480f275d7a96b0e3ddb", "67f2e1bc025833b0fa31f47195198997"};
        WalkerTestSpec spec7 = new WalkerTestSpec(
                "-T VariantFiltration -X OnOffGenotypeRatio:threshold=0.9 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5OnOff));
        executeTest("testOnOffGenotypeFilter", spec7);

        String[] md5Clusters = {"f5c7c5da0198c4aaafdfcfb3c356eedc", "8fa6b6ffc93ee7fb8d6b52a7fb7815ef"};
        WalkerTestSpec spec8 = new WalkerTestSpec(
                "-T VariantFiltration -X ClusteredSnps:window=10,snps=3 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Clusters));
        executeTest("testClusteredSnpsFilter", spec8);

        String[] md5Indels = {"30bf4c764f6dfd006d919ecaceee0166", "8e0e915a1cb63d7049e0671ed00101fe"};
        WalkerTestSpec spec9 = new WalkerTestSpec(
                "-T VariantFiltration -X IndelArtifact -B indels,PointIndel,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.indels -B cleaned,CleanedOutSNP,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.realigner_badsnps -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Indels));
        executeTest("testIndelArtifactFilter", spec9);






    }
}