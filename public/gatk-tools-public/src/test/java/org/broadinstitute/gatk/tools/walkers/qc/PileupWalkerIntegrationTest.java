/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.qc;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PileupWalkerIntegrationTest extends WalkerTest {

    @Test
    public void testGnarleyFHSPileup() {
        String gatk_args = "-T Pileup -I " + validationDataLocation + "FHS_Pileup_Test.bam "
                 + "-R " + hg18Reference
                 +  " -L chr15:46,347,148 -o %s";
        String expected_md5 = "526c93b0fa660d6b953b57103e59fe98";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(expected_md5));
        executeTest("Testing the standard (no-indel) pileup on three merged FHS pools with 27 deletions in 969 bases", spec);
    }



    private final static String SingleReadAligningOffChromosome1MD5 = "4a45fe1f85aaa8c4158782f2b6dee2bd";
    @Test
    public void testSingleReadAligningOffChromosome1() {
        String gatk_args = "-T Pileup "
                + " -I " + privateTestDir + "readOffb37contig1.bam"
                + " -R " + b37KGReference
                + " -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(SingleReadAligningOffChromosome1MD5));
        executeTest("Testing single read spanning off chromosome 1", spec);
    }

    @Test
    public void testSingleReadAligningOffChromosome1NoIndex() {
        String gatk_args = "-T Pileup "
                + " -I " + privateTestDir + "readOffb37contig1.noIndex.bam"
                + " -R " + b37KGReference
                + " -U ALLOW_UNINDEXED_BAM"
                + " -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(SingleReadAligningOffChromosome1MD5));
        executeTest("Testing single read spanning off chromosome 1 unindexed", spec);
    }

    /************************/

    //testing speedup to GATKBAMIndex


    @DataProvider(name="GATKBAMIndexTest")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();
        tests.add(new Object[]{"-L 20:1-76,050","8702701350de11a6d28204acefdc4775"});
        tests.add(new Object[]{"-L 20:10,000,000-10,001,100","818cf5a8229efe6f89fc1cd8145ccbe3"});
        tests.add(new Object[]{"-L 20:62,954,114-63,025,520","22471ea4a12e5139aef62bf8ff2a5b63"});
        tests.add(new Object[]{"-L 20:1-76,050 -L 20:20,000,000-20,000,100 -L 20:40,000,000-40,000,100 -L 20:30,000,000-30,000,100 -L 20:50,000,000-50,000,100 -L 20:62,954,114-63,025,520 ","08d899ed7c5a76ef3947bf67338acda1"});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GATKBAMIndexTest")
    public void testGATKBAMIndexSpeedup(final String intervals, final String md5){
        final String gatkArgs="-T Pileup -I " + validationDataLocation + "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam "
                + "-R " + hg19Reference + " -o %s ";

        WalkerTestSpec spec = new WalkerTestSpec(gatkArgs + intervals, 1, Arrays.asList(md5));
        executeTest("Testing with intervals="+intervals, spec);
    }


    /***********************/

    // testing hidden option -outputInsertLength
    private final static String SingleReadAligningOffChromosome1withInsertLengthMD5 = "279e2ec8832e540f47a6e2bdf4cef5ea";
    @Test
    public void testSingleReadAligningOffChromosome1withInsertLength() {
        String gatk_args = "-T Pileup "
                + " -I " + privateTestDir + "readOffb37contig1.bam"
                + " -R " + b37KGReference
                + " -outputInsertLength"
                + " -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(SingleReadAligningOffChromosome1withInsertLengthMD5));
        executeTest("Testing single read spanning off chromosome 1 (with insert length)", spec);
    }

    @Test
    public void testGnarleyFHSPileupwithInsertLength() {
        String gatk_args = "-T Pileup -I " + validationDataLocation + "FHS_Pileup_Test.bam "
                + "-R " + hg18Reference
                + " -outputInsertLength"
                +  " -L chr15:46,347,148 -o %s";
        String expected_md5 = "53ced173768f3d4d90b8a8206e72eae5";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 1, Arrays.asList(expected_md5));
        executeTest("Testing the standard (no-indel) pileup on three merged FHS pools with 27 deletions in 969 bases (with insert length)", spec);
    }

}
