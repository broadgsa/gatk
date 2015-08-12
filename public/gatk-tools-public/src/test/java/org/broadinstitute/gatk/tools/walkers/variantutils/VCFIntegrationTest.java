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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.AbstractIndex;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.interval.IntervalTreeIndex;
import htsjdk.tribble.index.linear.LinearIndex;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class VCFIntegrationTest extends WalkerTest {

    @Test(enabled = true)
    public void testReadingAndWritingWitHNoChanges() {

        String md5ofInputVCF = "3dc9ac85f2c0541df9bc57b4d81f480b";
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
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("b8b18a9a015cba1f3a44af532bf45338"));
        executeTest("Test reading and writing breakpoint VCF", spec1);
    }

    @Test(enabled = true)
    public void testReadingLowerCaseBases() {
        String testVCF = privateTestDir + "lowercaseBases.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("89daaa81f64e96cf4d9aa0abf9be0b76"));
        executeTest("Test reading VCF with lower-case bases", spec1);
    }

    @Test(enabled = true)
    public void testReadingAndWriting1000GSVs() {
        String testVCF = privateTestDir + "1000G_SVs.chr1.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("5d91e685c760f7e869cb06596d741116"));
        executeTest("Test reading and writing 1000G Phase I SVs", spec1);
    }

    @Test
    public void testReadingAndWritingSamtools() {
        String testVCF = privateTestDir + "samtools.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("122340e3dc333d2b4b79c5c0c443a3fe"));
        executeTest("Test reading and writing samtools vcf", spec1);
    }

    @Test
    public void testWritingSamtoolsWExBCFExample() {
        String testVCF = privateTestDir + "ex2.vcf";
        String baseCommand = "-R " + b36KGReference + " --no_cmdline_in_header -o %s ";
        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("db565efb14b2fe5f00a11762751d2476"));
        executeTest("Test writing samtools WEx BCF example", spec1);
    }

    @Test(enabled = true)
    public void testReadingSamtoolsWExBCFExample() {
        String testVCF = privateTestDir + "ex2.bcf";
        String baseCommand = "-R " + b36KGReference + " --no_cmdline_in_header -o %s ";
        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("0ca1a078d4801886ef4abac327df7104"));
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
        runVCFWithoutHeaders("-U LENIENT_VCF_PROCESSING", "a8f4be8ad9820286ea13a28a675133f1", null, true);
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
    public void testVCFIndexCreation(VCFIndexCreatorTest testSpec) throws NoSuchFieldException, IllegalAccessException, IOException {

        final String logFileName = new String("testVCFIndexCreation.log");
        final String chromosome = "20";
        final String commandLine = " -T SelectVariants" +
                " -R " + b37KGReference +
                " --no_cmdline_in_header" +
                " -L " + chromosome +
                " -V " + b37_NA12878_OMNI +
                " --variant_index_type " + testSpec.type +
                " --variant_index_parameter " + testSpec.parameter +
                " -log " + logFileName +
                " -o %s";
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        final String name = "testVCFIndexCreation: " + testSpec.toString();

        // execute that test and check if the actual and expected indices are the same
        executeTestAndCheckIndices(name, chromosome, testSpec, spec);

        // check the log for the warning message
        File file = new File(logFileName);
        Assert.assertTrue(FileUtils.readFileToString(file).contains(GATKVCFUtils.DEPRECATED_INDEX_ARGS_MSG));
    }

    @Test
    public void testVCFIndexCreationNoArgs() throws NoSuchFieldException, IllegalAccessException {

        final String chromosome = "20";
        final String commandLine = " -T SelectVariants" +
                " -R " + b37KGReference +
                " --no_cmdline_in_header" +
                " -L " + chromosome +
                " -V " + b37_NA12878_OMNI +
                " -o %s";
        final String name = "testVCFIndexCreationNoArgs";
        VCFIndexCreatorTest testSpec = new VCFIndexCreatorTest(GATKVCFUtils.DEFAULT_INDEX_TYPE, GATKVCFUtils.DEFAULT_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine, 1, Arrays.asList(""));
        spec.disableShadowBCF();

        // execute that test and check if the actual and expected indices are the same
        executeTestAndCheckIndices(name, chromosome, testSpec, spec);
    }

    @Test
    public void testGVCFIndexCreation() throws NoSuchFieldException, IllegalAccessException {

        final String chromosome = "20";
        final String commandLine = " -T SelectVariants" +
                " -R " + b37KGReference +
                " --no_cmdline_in_header" +
                " -L " + chromosome +
                " -V " + b37_NA12878_OMNI +
                " -o %s";
        final String name = "testGVCFIndexCreation";
        VCFIndexCreatorTest testSpec = new VCFIndexCreatorTest(GATKVCFUtils.DEFAULT_GVCF_INDEX_TYPE, GATKVCFUtils.DEFAULT_GVCF_INDEX_PARAMETER);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine, Arrays.asList(GATKVCFUtils.GVCF_EXT), Arrays.asList(""));
        spec.disableShadowBCF();

        // execute that test and check if the actual and expected indices are the same
        executeTestAndCheckIndices(name, chromosome, testSpec, spec);
    }

    //
    //
    // Block-Compressed Tabix Index Tests
    //
    //

    private class BlockCompressedIndexCreatorTest extends TestDataProvider {
        private final String extension;

        private BlockCompressedIndexCreatorTest(String extension) {
            super(BlockCompressedIndexCreatorTest.class);

            this.extension = extension;
        }

        public String toString() {
            return String.format("File extension %s", extension);
        }
    }

    @DataProvider(name = "BlockCompressedIndexDataProvider")
    public Object[][] blockCompressedIndexCreatorData() {
        for (final String extension : AbstractFeatureReader.BLOCK_COMPRESSED_EXTENSIONS)
            new BlockCompressedIndexCreatorTest(".vcf" + extension);

        return TestDataProvider.getTests(BlockCompressedIndexCreatorTest.class);
    }

    @Test(dataProvider = "BlockCompressedIndexDataProvider")
    public void testBlockCompressedIndexCreation(BlockCompressedIndexCreatorTest testSpec) throws NoSuchFieldException, IllegalAccessException {

        final String commandLine = " -T SelectVariants" +
                " -R " + b37KGReference +
                " --no_cmdline_in_header" +
                " -L 20" +
                " -V " + b37_NA12878_OMNI;
        final String name = "testBlockCompressedIndexCreation: " + testSpec.toString();

        File outVCF = createTempFile("testBlockCompressedIndexCreation", testSpec.extension);
        final WalkerTestSpec spec = new WalkerTestSpec(commandLine, 1, Arrays.asList(""));
        spec.disableShadowBCF();
        spec.setOutputFileLocation(outVCF);

        executeTest(name, spec);

        File outTribbleIdx = new File(outVCF.getAbsolutePath() + Tribble.STANDARD_INDEX_EXTENSION);
        Assert.assertFalse(outTribbleIdx.exists(), "testBlockCompressedIndexCreation: Want Tabix index but Tribble index exists: " + outTribbleIdx);

        File outTabixIdx = new File(outVCF.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
        final Index actualIndex = IndexFactory.loadIndex(outTabixIdx.toString());
        Assert.assertTrue(actualIndex instanceof TabixIndex, "testBlockCompressedIndexCreation: Want Tabix index but index is not Tabix: " + outTabixIdx);
    }

    //
    //
    // Block-Compressed Input Tests
    //
    //

    private class BlockCompressedInputTest extends TestDataProvider {
        private final String extension;

        private BlockCompressedInputTest(String extension) {
            super(BlockCompressedInputTest.class);

            this.extension = extension;
        }

        public String toString() {
            return String.format("File extension %s", extension);
        }
    }

    @DataProvider(name = "BlockCompressedInputDataProvider")
    public Object[][] blockCompressedInputData() {
        for (final String extension : AbstractFeatureReader.BLOCK_COMPRESSED_EXTENSIONS)
            new BlockCompressedInputTest(".vcf" + extension);

        return TestDataProvider.getTests(BlockCompressedInputTest.class);
    }

    @Test(dataProvider = "BlockCompressedInputDataProvider")
    public void testBlockCompressedInput(BlockCompressedInputTest testSpec) {

        File inputFile = new File(BaseTest.privateTestDir, "block_compressed_input_test" + testSpec.extension);
        final String commandLine = " -T SelectVariants" +
                " -R " + b37KGReference +
                " --no_cmdline_in_header" +
                " -V " + inputFile +
                " -o %s ";
        final String name = "testBlockCompressedInput: " + testSpec.toString();

        final WalkerTestSpec spec = new WalkerTestSpec(commandLine, 1, Arrays.asList("ce9c0bf31ee9452ac4a12a59d5814545"));

        executeTest(name, spec);
    }

    private void executeTestAndCheckIndices(final String name, final String chr, final VCFIndexCreatorTest testSpec, final WalkerTestSpec walkerTestSpec)
            throws NoSuchFieldException, IllegalAccessException {

        File outVCF = executeTest(name, walkerTestSpec).first.get(0);
        File outIdx = new File(outVCF.getAbsolutePath() + Tribble.STANDARD_INDEX_EXTENSION);

        final Index actualIndex = IndexFactory.loadIndex(outIdx.getAbsolutePath());
        final Index expectedIndex = testSpec.getIndex(outVCF);

        if (testSpec.type.equals("LINEAR"))
            Assert.assertTrue(actualIndex instanceof LinearIndex, "Index is not a LinearIndex");
        else if (testSpec.type.equals("INTERVAL"))
            Assert.assertTrue(actualIndex instanceof IntervalTreeIndex, "Index is not a IntervalTreeIndex");
        // dynamic indices ultimately resolve to one of LinearIndex or IntervalTreeIndex

        Assert.assertTrue(GATKVCFUtils.equivalentAbstractIndices((AbstractIndex) actualIndex, (AbstractIndex) expectedIndex), "Indices are not equivalent");

        if (actualIndex instanceof LinearIndex && expectedIndex instanceof LinearIndex) {
            Assert.assertTrue(GATKVCFUtils.equivalentLinearIndices((LinearIndex) actualIndex, (LinearIndex) expectedIndex, chr), "Linear indices are not equivalent");
        }
        else if (actualIndex instanceof IntervalTreeIndex && expectedIndex instanceof IntervalTreeIndex) {
            Assert.assertTrue(GATKVCFUtils.equivalentIntervalIndices((IntervalTreeIndex) actualIndex, (IntervalTreeIndex) expectedIndex, chr), "Interval indices are not equivalent");
        }
        else {
            Assert.fail("Indices are not of the same type");
        }
    }
}
