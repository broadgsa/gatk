package org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator;


import java.util.Arrays;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

public class
        GenomicAnnotatorIntegrationTest extends WalkerTest {
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


        String[] md5WithDashSArg = {"9583d7060bc3de73b392e7435c31946b"};
        WalkerTestSpec specWithSArg = new WalkerTestSpec(
                "-T GenomicAnnotator -R " + b36KGReference +
                " -B:variant,vcf /humgen/gsa-hpprojects/GATK/data/Annotations/examples/CEU_hapmap_nogt_23_subset.vcf" +
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
}
