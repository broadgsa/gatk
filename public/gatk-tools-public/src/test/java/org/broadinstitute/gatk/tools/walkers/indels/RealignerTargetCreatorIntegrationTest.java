/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class RealignerTargetCreatorIntegrationTest extends WalkerTest {

    @DataProvider(name = "intervals1")
    public Object[][] intervals1() {
        String arguments = "-T RealignerTargetCreator -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam --mismatchFraction 0.15 -L 1:10,000,000-10,050,000";
        return new Object[][]{
                {"test standard nt=1", arguments},
                {"test standard nt=4", "-nt 4 " + arguments}
        };
    }

    @DataProvider(name = "intervals2")
    public Object[][] intervals2() {
        String arguments = "-T RealignerTargetCreator --known " + b36dbSNP129 + " -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,200,000";
        return new Object[][]{
                {"test with dbsnp nt=1", arguments},
                {"test with dbsnp nt=4", "-nt 4 " + arguments}
        };
    }

    @Test(dataProvider = "intervals1")
    public void testIntervals1(String testName, String arguments) {
        String md5 = "3f0b63a393104d0c4158c7d1538153b8";

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(arguments + " -o %s", 1, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test(dataProvider = "intervals2")
    public void testIntervals2(String testName, String arguments) {
        String md5 = "d073237694175c75d37bd4f40b8c64db";

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(arguments + " -o %s", 1, Arrays.asList(md5));
        executeTest(testName, spec);
    }

    @Test
    public void testKnownsOnly() {
        WalkerTest.WalkerTestSpec spec3 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -R " + b36KGReference + " --known " + privateTestDir + "NA12878.chr1_10mb_11mb.slx.indels.vcf4 -L " + privateTestDir + "NA12878.chr1_10mb_11mb.slx.indels.vcf4 -o %s",
                 1,
                 Arrays.asList("5206cee6c01b299417bf2feeb8b3dc96"));
        executeTest("test rods only", spec3);
    }

    @Test()
    public void testBadCigarStringDoesNotFail() {
        // Just making sure the test runs without an error, don't care about the MD5 value
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T RealignerTargetCreator -R " + b37KGReference + " -I " + privateTestDir + "Realigner.error.bam -L 19:5787200-5787300 -o %s",
                1,
                Arrays.asList(""));
        executeTest("test bad cigar string string does not fail", spec);
    }

    @Test(dataProvider = "intervals1")
    public void testTargetListAgainstIntervalList(String testName, String arguments) throws IOException {
        final List<String> md5 = Collections.emptyList();
        final File targetListFile = createTempFile("RTCTest", ".targets");
        final File intervalListFile = createTempFile("RTCTest", ".interval_list");

        WalkerTest.WalkerTestSpec targetListSpec = new WalkerTest.WalkerTestSpec(arguments, 1, md5);
        WalkerTest.WalkerTestSpec intervalListSpec = new WalkerTest.WalkerTestSpec(arguments, 1, md5);

        targetListSpec.setOutputFileLocation(targetListFile);
        intervalListSpec.setOutputFileLocation(intervalListFile);

        executeTest(testName + " (compare target-list and interval-list output)", targetListSpec);
        executeTest(testName + " (compare target-list and interval-list output)", intervalListSpec);

        final ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(BaseTest.hg19Reference));
        final GenomeLocParser hg19GenomeLocParser = new GenomeLocParser(seq);
        final List<GenomeLoc> targetList = IntervalUtils.intervalFileToList(hg19GenomeLocParser,
                targetListFile.getAbsolutePath());
        final List<Interval> targetListResult = new ArrayList<>();
        for ( GenomeLoc target : targetList ) {
            targetListResult.add(new Interval(target.getContig(), target.getStart(), target.getStop()));
        }

        final List<Interval> intervalListResult = IntervalList.fromFile(intervalListFile).getIntervals();

        Assert.assertFalse(targetListResult.isEmpty());
        Assert.assertFalse(intervalListResult.isEmpty());
        Assert.assertEquals(targetListResult, intervalListResult);
    }
}
