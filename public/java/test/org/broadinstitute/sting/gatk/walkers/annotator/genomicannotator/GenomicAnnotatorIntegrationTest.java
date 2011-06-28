package org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator;


import java.util.Arrays;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

public class GenomicAnnotatorIntegrationTest extends WalkerTest {
    String testFileWithIndels = validationDataLocation + "/GenomicAnnotatorValidation/1KGBroadWEx.cleaned.indels.vcf";
    String testFileWithSNPsAndIndels = validationDataLocation + "/GenomicAnnotatorValidation/1KGBroadWEx.variants.vcf";

    @Test
    public void testGenomicAnnotatorOnDbSNP() {

        /*
        TODO put this test back in once it gets faster.
        String[] md5 = {"d19d6d1eb52fb09e7493653dc645d92a"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T GenomicAnnotator -R " + b36KGReference + " " +
                "-B:variant,vcf /humgen/gsa-hpprojects/GATK/data/Annotations/examples/CEU_hapmap_nogt_23_subset.vcf " +
                "-B:dbsnp,AnnotatorInputTable /humgen/gsa-hpprojects/GATK/data/Annotations/dbsnp/b130/snp130-b36-only-the-SNPs.txt " +
                "-m " + //generate many records from one input record if necessary
                "-o %s " +
                "-BTI variant",
                 1,
                 Arrays.asList(md5));
        executeTest("test with dbSNP", spec);
        */


        String[] md5WithDashSArg = {"3d3b61a83c1189108eabb2df04218099"};
        WalkerTestSpec specWithSArg = new WalkerTestSpec(
                "-T GenomicAnnotator -R " + b36KGReference +
                " -B:variant,vcf3 /humgen/gsa-hpprojects/GATK/data/Annotations/examples/CEU_hapmap_nogt_23_subset.vcf" +
                " -B:dbsnp,AnnotatorInputTable /humgen/gsa-hpprojects/GATK/data/Annotations/dbsnp/b130/snp130-b36-only-the-SNPs.txt" +
                " -m" + //generate many records from one input record if necessary
                " -o %s" +
                " -BTI variant" +
                " -s dbsnp.name,dbsnp.refUCSC,dbsnp.strand,dbsnp.observed,dbsnp.avHet" +
                " -NO_HEADER",
                 1,
                 Arrays.asList(md5WithDashSArg));
        executeTest("test with dbSNP and -s arg", specWithSArg);

    }

    @Test
    public void testGenomicAnnotatorOnIndels() {
        WalkerTestSpec testOnIndels = new WalkerTestSpec(
                buildCommandLine(
                        "-T GenomicAnnotator",
                        "-R " + b37KGReference,
                        "-L 22:10000000-20000000",
                        "-B:refseq,AnnotatorInputTable " + b37Refseq,
                        "-B:variant,VCF " + testFileWithIndels,
                        "-NO_HEADER",
                        "-o %s"
                ),
                1,
                Arrays.asList("caa562160733aa638e1ba413ede209ae")
        );
        executeTest("testGenomicAnnotatorOnIndels", testOnIndels);
    }

    @Test
    public void testGenomicAnnotatorOnSNPsAndIndels() {
        WalkerTestSpec testOnSNPsAndIndels = new WalkerTestSpec(
                buildCommandLine(
                        "-T GenomicAnnotator",
                        "-R " + b37KGReference,
                        "-L 22:10000000-20000000",
                        "-B:refseq,AnnotatorInputTable " + b37Refseq,
                        "-B:variant,VCF " + testFileWithSNPsAndIndels,
                        "-NO_HEADER",
                        "-o %s"
                ),
                1,
                Arrays.asList("a4cf76f08fa90284b6988a464b6e0c17")
        );
        executeTest("testGenomicAnnotatorOnSNPsAndIndels", testOnSNPsAndIndels);
    }
}
