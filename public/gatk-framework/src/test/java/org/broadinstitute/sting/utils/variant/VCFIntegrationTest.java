/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.variant;

import org.broad.tribble.index.AbstractIndex;
import org.broad.tribble.index.ChrIndex;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.index.interval.IntervalTreeIndex;
import org.broad.tribble.index.linear.LinearIndex;
import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

public class VCFIntegrationTest extends WalkerTest {

    @Test(enabled = true)
    public void testReadingAndWritingWitHNoChanges() {

        String md5ofInputVCF = "d991abe6c6a7a778a60a667717903be0";
        String testVCF = privateTestDir + "vcf4.1.example.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T VariantAnnotator --variant " + testVCF + " -L " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList(md5ofInputVCF));
        List<File> result = executeTest("Test Variant Annotator with no changes", spec1).getFirst();

        String test2 = baseCommand + "-T VariantsToVCF --variant " + result.get(0).getAbsolutePath();
        WalkerTestSpec spec2 = new WalkerTestSpec(test2, 1, Arrays.asList(md5ofInputVCF));
        executeTest("Test Variants To VCF from new output", spec2);
    }

    @Test(enabled = true)
    public void testReadingAndWritingBreakpointAlleles() {
        String testVCF = privateTestDir + "breakpoint-example.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("13329ba7360a8beb3afc02569e5a20c4"));
        executeTest("Test reading and writing breakpoint VCF", spec1);
    }

    @Test(enabled = true)
    public void testReadingLowerCaseBases() {
        String testVCF = privateTestDir + "lowercaseBases.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("e0e308a25e56bde1c664139bb44ed19d"));
        executeTest("Test reading VCF with lower-case bases", spec1);
    }

    @Test(enabled = true)
    public void testReadingAndWriting1000GSVs() {
        String testVCF = privateTestDir + "1000G_SVs.chr1.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("bdab26dd7648a806dbab01f64db2bdab"));
        executeTest("Test reading and writing 1000G Phase I SVs", spec1);
    }

    @Test
    public void testReadingAndWritingSamtools() {
        String testVCF = privateTestDir + "samtools.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("38697c195e7abf18d95dcc16c8e6d284"));
        executeTest("Test reading and writing samtools vcf", spec1);
    }

    @Test
    public void testWritingSamtoolsWExBCFExample() {
        String testVCF = privateTestDir + "ex2.vcf";
        String baseCommand = "-R " + b36KGReference + " --no_cmdline_in_header -o %s ";
        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("e8f721ce81e4fdadba13c5291027057f"));
        executeTest("Test writing samtools WEx BCF example", spec1);
    }

    @Test(enabled = true)
    public void testReadingSamtoolsWExBCFExample() {
        String testVCF = privateTestDir + "ex2.bcf";
        String baseCommand = "-R " + b36KGReference + " --no_cmdline_in_header -o %s ";
        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("0439e2b4ccc63bb4ba7c283cd9ab1b25"));
        executeTest("Test reading samtools WEx BCF example", spec1);
    }

    //
    //
    // Tests to ensure that -U LENIENT_VCF_PROCESS
    //
    //

    @Test
    public void testFailingOnVCFWithoutHeaders() {
        runVCFWithoutHeaders("", "", IllegalStateException.class, false);
    }

    @Test
    public void testPassingOnVCFWithoutHeadersWithLenientProcessing() {
        runVCFWithoutHeaders("-U LENIENT_VCF_PROCESSING", "6de8cb7457154dd355aa55befb943f88", null, true);
    }

    private void runVCFWithoutHeaders(final String moreArgs, final String expectedMD5, final Class expectedException, final boolean disableBCF) {
        final String testVCF = privateTestDir + "vcfexample2.noHeader.vcf";
        final String baseCommand = "-R " + b37KGReference
                + " --no_cmdline_in_header -o %s "
                + "-T VariantsToVCF -V " + testVCF + " " + moreArgs;
        WalkerTestSpec spec1 = expectedException != null
                ? new WalkerTestSpec(baseCommand, 1, expectedException)
                : new WalkerTestSpec(baseCommand, 1, Arrays.asList(expectedMD5));
        if ( disableBCF )
            spec1.disableShadowBCF();
        executeTest("Test reading VCF without header lines with additional args " + moreArgs, spec1);
    }

    //
    //
    // IndexCreator tests
    //
    //

    private class VCFIndexCreatorTest extends TestDataProvider {
        private final GATKVCFIndexType type;
        private final int parameter;

        private VCFIndexCreatorTest(GATKVCFIndexType type, int parameter) {
            super(VCFIndexCreatorTest.class);

            this.type = type;
            this.parameter = parameter;
        }

        public String toString() {
            return String.format("Index Type %s, Index Parameter %s", type, parameter);
        }

        public Index getIndex(final File vcfFile) {
            switch (type) {
                case DYNAMIC_SEEK : return IndexFactory.createDynamicIndex(vcfFile, new VCFCodec(), IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
                case DYNAMIC_SIZE : return IndexFactory.createDynamicIndex(vcfFile, new VCFCodec(), IndexFactory.IndexBalanceApproach.FOR_SIZE);
                case LINEAR : return IndexFactory.createLinearIndex(vcfFile, new VCFCodec(), parameter);
                case INTERVAL : return IndexFactory.createIntervalIndex(vcfFile, new VCFCodec(), parameter);
                default : throw new TestException("Invalid index type");
            }
        }
    }

    @DataProvider(name = "IndexDataProvider")
    public Object[][] indexCreatorData() {
        new VCFIndexCreatorTest(GATKVCFIndexType.DYNAMIC_SEEK, 0);
        new VCFIndexCreatorTest(GATKVCFIndexType.DYNAMIC_SIZE, 0);
        new VCFIndexCreatorTest(GATKVCFIndexType.LINEAR, 100);
        new VCFIndexCreatorTest(GATKVCFIndexType.LINEAR, 10000);
        new VCFIndexCreatorTest(GATKVCFIndexType.INTERVAL, 20);
        new VCFIndexCreatorTest(GATKVCFIndexType.INTERVAL, 2000);

        return TestDataProvider.getTests(VCFIndexCreatorTest.class);
    }

    @Test(dataProvider = "IndexDataProvider")
    public void testVCFIndexCreation(VCFIndexCreatorTest testSpec) throws NoSuchFieldException, IllegalAccessException {

        final String commandLine = " -T SelectVariants" +
                " -R " + b37KGReference +
                " --no_cmdline_in_header" +
                " -L 20" +
                " -V " + b37_NA12878_OMNI +
                " --variant_index_type " + testSpec.type +
                " --variant_index_parameter " + testSpec.parameter +
                " -o %s ";
        final String name = "testVCFIndexCreation: " + testSpec.toString();

        final WalkerTestSpec spec = new WalkerTestSpec(commandLine, 1, Arrays.asList(""));
        spec.disableShadowBCF();

        File outVCF = executeTest(name, spec).first.get(0);
        File outIdx = new File(outVCF.getAbsolutePath() + ".idx");

        final Index actualIndex = IndexFactory.loadIndex(outIdx.getAbsolutePath());
        final Index expectedIndex = testSpec.getIndex(outVCF);

        if (testSpec.type.equals("LINEAR"))
            Assert.assertTrue(actualIndex instanceof LinearIndex, "Index is not a LinearIndex");
        else if (testSpec.type.equals("INTERVAL"))
            Assert.assertTrue(actualIndex instanceof IntervalTreeIndex, "Index is not a IntervalTreeIndex");
        // dynamic indices ultimately resolve to one of LinearIndex or IntervalTreeIndex

        Assert.assertTrue(equivalentAbstractIndices((AbstractIndex)actualIndex, (AbstractIndex)expectedIndex), "Indices are not equivalent");

        if (actualIndex instanceof LinearIndex && expectedIndex instanceof LinearIndex) {
            Assert.assertTrue(equivalentLinearIndices((LinearIndex)actualIndex, (LinearIndex)expectedIndex, "20"), "Linear indices are not equivalent");
        }
        else if (actualIndex instanceof IntervalTreeIndex && expectedIndex instanceof IntervalTreeIndex) {
            Assert.assertTrue(equivalentIntervalIndices((IntervalTreeIndex)actualIndex, (IntervalTreeIndex)expectedIndex, "20"), "Interval indices are not equivalent");
        }
        else {
            Assert.fail("Indices are not of the same type");
        }
    }

    private static boolean equivalentAbstractIndices(AbstractIndex thisIndex, AbstractIndex otherIndex){
        return thisIndex.getVersion() == otherIndex.getVersion() &&
                thisIndex.getIndexedFile().equals(otherIndex.getIndexedFile()) &&
                thisIndex.getIndexedFileSize() == otherIndex.getIndexedFileSize() &&
                thisIndex.getIndexedFileMD5().equals(otherIndex.getIndexedFileMD5()) &&
                thisIndex.getFlags() == otherIndex.getFlags();
     }

    private static boolean equivalentLinearIndices(LinearIndex thisIndex, LinearIndex otherIndex, String chr) throws NoSuchFieldException, IllegalAccessException {
        org.broad.tribble.index.linear.LinearIndex.ChrIndex thisChr = (org.broad.tribble.index.linear.LinearIndex.ChrIndex)getChrIndex(thisIndex, chr);
        org.broad.tribble.index.linear.LinearIndex.ChrIndex otherChr = (org.broad.tribble.index.linear.LinearIndex.ChrIndex)getChrIndex(otherIndex, chr);

        return  thisChr.getName().equals(otherChr.getName()) &&
                //thisChr.getTotalSize() == otherChr.getTotalSize() &&      TODO: why does this differ?
                thisChr.getNFeatures() == otherChr.getNFeatures() &&
                thisChr.getNBlocks() == otherChr.getNBlocks();
    }

    private static boolean equivalentIntervalIndices(IntervalTreeIndex thisIndex, IntervalTreeIndex otherIndex, String chr) throws NoSuchFieldException, IllegalAccessException {
        org.broad.tribble.index.interval.IntervalTreeIndex.ChrIndex thisChr = (org.broad.tribble.index.interval.IntervalTreeIndex.ChrIndex)getChrIndex(thisIndex, chr);
        org.broad.tribble.index.interval.IntervalTreeIndex.ChrIndex otherChr = (org.broad.tribble.index.interval.IntervalTreeIndex.ChrIndex)getChrIndex(otherIndex, chr);

        // TODO: compare trees?
        return thisChr.getName().equals(otherChr.getName());
    }

    private static ChrIndex getChrIndex(AbstractIndex index, String chr) throws NoSuchFieldException, IllegalAccessException {
        Field f = AbstractIndex.class.getDeclaredField("chrIndices");
        f.setAccessible(true);
        LinkedHashMap<String, ChrIndex> chrIndices = (LinkedHashMap<String, ChrIndex>) f.get(index);
        return chrIndices.get(chr);
    }
}
