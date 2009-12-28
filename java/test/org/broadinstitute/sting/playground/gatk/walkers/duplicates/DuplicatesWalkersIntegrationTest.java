package org.broadinstitute.sting.playground.gatk.walkers.duplicates;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class DuplicatesWalkersIntegrationTest extends WalkerTest {
    public void testClipper(String name, String args, String md5) {
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

    @Test public void testChr110Mb() { testClipper("testChr1-10mb", "-L chr1:1-10,000,000 --quiet", ""); }
    @Test public void testIntervalVerbose() { testClipper("testIntervalVerbose", "-L chr1:6,527,154-6,528,292", ""); }
}