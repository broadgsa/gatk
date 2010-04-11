package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.ArrayList;

public class UnifiedGenotyperPerformanceTest extends WalkerTest {

    @Test
    public void testUnifiedGenotyperWholeGenome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " -T UnifiedGenotyper" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-50,000,000" +
                        " -varout /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testUnifiedGenotyperWholeGenome", spec);
    }

    @Test
    public void testUnifiedGenotyperWholeExome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " -T UnifiedGenotyper" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L " + evaluationDataLocation + "whole_exome_agilent_designed_120.targets.chr1.interval_list" +
                        " -varout /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testUnifiedGenotyperWholeExome", spec);
    }
}