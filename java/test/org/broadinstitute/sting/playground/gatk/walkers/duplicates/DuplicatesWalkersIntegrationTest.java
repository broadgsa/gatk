package org.broadinstitute.sting.playground.gatk.walkers.duplicates;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class DuplicatesWalkersIntegrationTest extends WalkerTest {
    public void testCounter(String name, String args, String md5) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CountDuplicates " +
                        "-R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta " +
                        "-I /humgen/gsa-hpprojects/GATK/data/Validation_Data/TCGA-06-0188.aligned.duplicates_marked.bam " +
                        "-o %s " + args,
                1, // just one output file
                Arrays.asList("tmp"),
                Arrays.asList(md5));
        List<File> result = executeTest(name, spec).getFirst();
    }

    @Test public void testChr110Mb() { testCounter("testChr1-10mb", "-L chr1:1-10,000,000 --quiet", "fa8bfdd0b62a13a543bae90f7c674db7"); }
    @Test public void testIntervalVerbose() { testCounter("testIntervalVerbose", "-L chr1:6,527,154-6,528,292", "1ebcc10b85af16805a54391721776657"); }

    public void testCombiner(String name, String args, String md51, String md52) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T CombineDuplicates " +
                        "-R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta " +
                        "-I /humgen/gsa-hpprojects/GATK/data/Validation_Data/TCGA-06-0188.aligned.duplicates_marked.bam " +
                        "-o %s --outputBAM %s " + args,
                2, // just one output file
                Arrays.asList("tmp", "bam"),
                Arrays.asList(md51, md52));
        List<File> result = executeTest(name, spec).getFirst();
    }

    @Test public void testIntervalCombine() { testCombiner("testIntervalCombine", "-L chr1:6,527,154-6,528,292 -maxQ 50", "d41d8cd98f00b204e9800998ecf8427e", "e2501d7e20564df9f7519fadef8cf283"); }
    @Test public void testIntervalCombineQ60() { testCombiner("testIntervalCombine", "-L chr1:6,527,154-6,528,292 -maxQ 60", "d41d8cd98f00b204e9800998ecf8427e", "b23f6436d230f57f969502ddd8d48c18"); }
}