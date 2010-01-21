package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 11, 2010
 * Time: 11:30:59 AM
 * To change this template use File | Settings | File Templates.
 */
public class SequenomToVCFIntegrationTest extends WalkerTest {
    @Test
    public void testSequenomToVCFWithoutSmartHardy() {
        String testMD5 = "a16ce402ce4492c8c0552073c0a8bdf5";
        String testPedFile = "/humgen/gsa-hpprojects/GATK/data/Validation_Data/Sequenom_Test_File.txt";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -L 1:1000000-2000000 "+
                          "-T PlinkToVCF -b36contig -ns 10 -sPed "+testPedFile+" -vcf %s";
        WalkerTest.WalkerTestSpec spec = new WalkerTestSpec(testArgs,1, Arrays.asList(testMD5));
        List<File> result = executeTest("TestSequenomToVCFNoSmartHardy",spec).getFirst();
    }
}
