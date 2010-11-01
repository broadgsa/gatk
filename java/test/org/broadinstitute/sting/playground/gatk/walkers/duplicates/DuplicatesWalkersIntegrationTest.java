package org.broadinstitute.sting.playground.gatk.walkers.duplicates;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class DuplicatesWalkersIntegrationTest extends WalkerTest {
    public void testCounter(String name, String args, String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CountDuplicates" +
                        " -R " + hg18Reference +
                        " -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/TCGA-06-0188.aligned.duplicates_marked.bam" +
                        " -o %s " + args,
                1, // just one output file
                Arrays.asList("tmp"),
                Arrays.asList(md5));
        List<File> result = executeTest(name, spec).getFirst();
    }

    @Test public void testChr110Mb() { testCounter("testChr1-10mb", "-L chr1:1-10,000,000 --quietLocus", "d3c329a634904d95c4b180d0d63eadfc"); }
    @Test public void testIntervalVerbose() { testCounter("testIntervalVerbose", "-L chr1:6,527,154-6,528,292", "5fbb930020df6ca7d0f724524fc43b3e"); }

    public void testCombiner(String name, String args, String md51, String md52) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineDuplicates" +
                        " -R " + hg18Reference +
                        " -I /humgen/gsa-hpprojects/GATK/data/Validation_Data/TCGA-06-0188.aligned.duplicates_marked.bam" +
                        " -o %s --outputBAM %s " + args,
                2, // just one output file
                Arrays.asList("tmp", "bam"),
                Arrays.asList(md51, md52));
        List<File> result = executeTest(name, spec).getFirst();
    }

    @Test public void testIntervalCombine() { testCombiner("testIntervalCombine", "-L chr1:6,527,154-6,528,292 -maxQ 50", "d41d8cd98f00b204e9800998ecf8427e", "4541f57820637039bc2f5a97bcaadfe4"); }
    @Test public void testIntervalCombineQ60() { testCombiner("testIntervalCombine", "-L chr1:6,527,154-6,528,292 -maxQ 60", "d41d8cd98f00b204e9800998ecf8427e", "8c0350c0a697e4083aab6ead3f404de4"); }
}
