package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {
        String[] md5DoC = {"d11547059193b90e2fd1eb47dac18999", "21c8e1f9dc65fdfb39347547f9b04011"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T VariantFiltration -X DepthOfCoverage:max=70 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5DoC));
        executeTest("testDoCFilter", spec1);

        String[] md5AlleleBalance = {"409f2b7fb193919d28c20be4dc5c4516", "a13e4ce6260bf9f33ca99dc808b8e6ad"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T VariantFiltration -X AlleleBalance:low=0.25,high=0.75 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5AlleleBalance));
        executeTest("testAlleleBalanceFilter", spec2);

        String[] md5Strand = {"e5806925867088edeee9f6276ab030f3", "0f7db0aad764268ee8fa3b857df8d87d"};
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T VariantFiltration -X FisherStrand:pvalue=0.0001 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Strand));
        executeTest("testStrandFilter", spec3);

        String[] md5Lod = {"2c93f56e07e320ac7a012bd4940b1dbc", "7e0c4f2b0fda85fd2891eee76c396a55"};
        WalkerTestSpec spec4 = new WalkerTestSpec(
                "-T VariantFiltration -X LodThreshold:lod=10 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Lod));
        executeTest("testLodFilter", spec4);

        String[] md5MQ0 = {"fde4be5eec12be3405863649170aaaac", "3203de335621851bccf596242b079e23"};
        WalkerTestSpec spec5 = new WalkerTestSpec(
                "-T VariantFiltration -X MappingQualityZero:max=70 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5MQ0));
        executeTest("testMappingQuality0Filter", spec5);

        String[] md5MQ = {"17c35a5ccd814991dbcd3be6c2f244cd", "ecc777feedea61f7b570d114c2ab89b1"};
        WalkerTestSpec spec6 = new WalkerTestSpec(
                "-T VariantFiltration -X MappingQuality:min=20 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5MQ));
        executeTest("testRMSMappingQualityFilter", spec6);

        String[] md5OnOff = {"fc649a7ebd4a7b7fa739fc4fcb7704e2", "67f2e1bc025833b0fa31f47195198997"};
        WalkerTestSpec spec7 = new WalkerTestSpec(
                "-T VariantFiltration -X OnOffGenotypeRatio:threshold=0.9 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5OnOff));
        executeTest("testOnOffGenotypeFilter", spec7);

        String[] md5Clusters = {"9b3aede689a25431ed2e1aac9dd73dc6", "8fa6b6ffc93ee7fb8d6b52a7fb7815ef"};
        WalkerTestSpec spec8 = new WalkerTestSpec(
                "-T VariantFiltration -X ClusteredSnps:window=10,snps=3 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Clusters));
        executeTest("testClusteredSnpsFilter", spec8);

        String[] md5Indels = {"c36cb0021f90b3509e483b57f4f41f6e", "8e0e915a1cb63d7049e0671ed00101fe"};
        WalkerTestSpec spec9 = new WalkerTestSpec(
                "-T VariantFiltration -X IndelArtifact -B indels,PointIndel,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.indels -B cleaned,CleanedOutSNP,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.realigner_badsnps -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -B variant,Variants,/humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chr1_10mb_11mb.slx.geli.calls -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod -vcf %s -included %s -sample NA12878",
                 2,
                 Arrays.asList(md5Indels));
        executeTest("testIndelArtifactFilter", spec9);






    }
}