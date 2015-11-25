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

package org.broadinstitute.gatk.engine;

import htsjdk.samtools.*;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.filters.MappingQualityUnavailableFilter;
import org.broadinstitute.gatk.engine.filters.DuplicateReadFilter;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.utils.variant.VCIterable;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;

/**
 *
 */
public class EngineFeaturesIntegrationTest extends WalkerTest {
    private void testBadRODBindingInput(String type, String name, Class c) {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintVariantsWalker -L 1:1 --variant:variant," + type + " "
                + b37dbSNP132 + " -R " + b37KGReference + " -o %s",
                1, c);
        executeTest(name, spec);
    }

    @Test() private void testBadRODBindingInputType1() {
        testBadRODBindingInput("beagle", "BEAGLE input to VCF expecting walker", UserException.BadArgumentValue.class);
    }

    @Test() private void testBadRODBindingInputType3() {
        testBadRODBindingInput("bed", "Bed input to VCF expecting walker", UserException.BadArgumentValue.class);
    }

    @Test() private void testBadRODBindingInputTypeUnknownType() {
        testBadRODBindingInput("bedXXX", "Unknown input to VCF expecting walker", UserException.UnknownTribbleType.class);
    }

    private void testMissingFile(String name, String missingBinding) {
        WalkerTestSpec spec = new WalkerTestSpec(missingBinding + " -R " + b37KGReference + " -o %s",
                1, UserException.CouldNotReadInputFile.class);
        executeTest(name, spec);
    }

    @Test() private void testMissingBAMnt1() {
        testMissingFile("missing BAM", "-T TestPrintReadsWalker -I missing.bam -nt 1");
    }
    @Test() private void testMissingBAMnt4() {
        testMissingFile("missing BAM", "-T TestPrintReadsWalker -I missing.bam -nt 4");
    }
    @Test() private void testMissingVCF() {
        testMissingFile("missing VCF", "-T TestPrintVariantsWalker -V missing.vcf");
    }
    @Test() private void testMissingInterval() {
        testMissingFile("missing interval", "-T TestPrintReadsWalker -L missing.interval_list -I " + b37GoodBAM);
    }


    // --------------------------------------------------------------------------------
    //
    // Test that our exceptions are coming back as we expect
    //
    // --------------------------------------------------------------------------------

    private class EngineErrorHandlingTestProvider extends TestDataProvider {
        final Class expectedException;
        final String args;
        final int iterationsToTest;

        public EngineErrorHandlingTestProvider(Class exceptedException, final String args) {
            super(EngineErrorHandlingTestProvider.class);
            this.expectedException = exceptedException;
            this.args = args;
            this.iterationsToTest = args.equals("") ? 1 : 10;
            setName(String.format("Engine error handling: expected %s with args %s", exceptedException, args));
        }
    }

    @DataProvider(name = "EngineErrorHandlingTestProvider")
    public Object[][] makeEngineErrorHandlingTestProvider() {
        for ( final FailMethod failMethod : FailMethod.values() ) {
            if ( failMethod == FailMethod.TREE_REDUCE )
                continue; // cannot reliably throw errors in TREE_REDUCE

            final String failArg = " -fail " + failMethod.name();
            for ( final String args : Arrays.asList("", " -nt 2", " -nct 2") ) {
                new EngineErrorHandlingTestProvider(NullPointerException.class, failArg + args);
                new EngineErrorHandlingTestProvider(UserException.class, failArg + args);
                new EngineErrorHandlingTestProvider(ReviewedGATKException.class, failArg + args);
            }
        }

        return EngineErrorHandlingTestProvider.getTests(EngineErrorHandlingTestProvider.class);
    }

    //
    // Loop over errors to throw, make sure they are the errors we get back from the engine, regardless of NT type
    //
    @Test(enabled = true, dataProvider = "EngineErrorHandlingTestProvider", timeOut = 60 * 1000 )
    public void testEngineErrorHandlingTestProvider(final EngineErrorHandlingTestProvider cfg) {
        for ( int i = 0; i < cfg.iterationsToTest; i++ ) {
            final String root = "-T TestErrorThrowingWalker -R " + exampleFASTA;
            final String args = root + cfg.args + " -E " + cfg.expectedException.getSimpleName();
            WalkerTestSpec spec = new WalkerTestSpec(args, 0, cfg.expectedException);

            executeTest(cfg.toString(), spec);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Test that read filters are being applied in the order we expect
    //
    // --------------------------------------------------------------------------------

    @ReadFilters({MappingQualityUnavailableFilter.class, DuplicateReadFilter.class})
    @DisabledReadFilters({DuplicateReadFilter.class})
    public static class DummyReadWalkerWithFilters extends ReadWalker<Integer, Integer> {
        @Output
        PrintStream out;

        @Override
        public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
            return 1;
        }

        @Override
        public Integer reduceInit() {
            return 0;
        }

        @Override
        public Integer reduce(Integer value, Integer sum) {
            return value + sum;
        }

        @Override
        public void onTraversalDone(Integer result) {
            out.println(result);
        }
    }

    @Test(enabled = true)
    public void testUserReadFilterAppliedBeforeWalker() {
        WalkerTestSpec spec = new WalkerTestSpec("-R " + b37KGReference + " -I " + privateTestDir + "allMAPQ255.bam"
                + " -T DummyReadWalkerWithFilters -o %s -L MT -rf ReassignMappingQuality",
                1, Arrays.asList("ecf27a776cdfc771defab1c5d19de9ab"));
        executeTest("testUserReadFilterAppliedBeforeWalker", spec);
    }

    @Test(enabled = true)
    public void testUserReadFilterDisabledAppliedBeforeWalker() {
        WalkerTestSpec spec = new WalkerTestSpec("-R " + b37KGReference + " -I " + privateTestDir + "allMAPQ255.bam"
                + " -T DummyReadWalkerWithFilters -o %s -L MT -drf DuplicateRead",
                1, Arrays.asList("897316929176464ebc9ad085f31e7284"));
        executeTest("testUserReadFilterDisabledAppliedBeforeWalker", spec);
    }

    @Test( enabled = true, expectedExceptions = RuntimeException.class )
    public void testUserReadFilterDisabledAppliedBeforeWalkerException() {
        WalkerTestSpec spec = new WalkerTestSpec("-R " + b37KGReference + " -I " + privateTestDir + "allMAPQ255.bam"
                + " -T DummyReadWalkerWithFilters -o %s -L MT -drf ReassignMappingQuality",
                1, Arrays.asList(""));
        executeTest("testUserReadFilterDisabledAppliedBeforeWalkerException", spec);
    }

    @Test
    public void testNegativeCompress() {
        testBadCompressArgument(-1);
    }

    @Test
    public void testTooBigCompress() {
        testBadCompressArgument(100);
    }

    private void testBadCompressArgument(final int compress) {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintReadsWalker -R " + b37KGReference + " -I " + privateTestDir + "NA12878.1_10mb_2_10mb.bam -o %s -compress " + compress,
                1, UserException.class);
        executeTest("badCompress " + compress, spec);
    }

    // --------------------------------------------------------------------------------
    //
    // Test that the VCF version key is what we expect
    //
    // --------------------------------------------------------------------------------
    @Test(enabled = true)
    public void testGATKVersionInVCF() throws Exception {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintVariantsWalker -R " + b37KGReference +
                " -V " + privateTestDir + "NA12878.WGS.b37.chr20.firstMB.vcf"
                + " -o %s -L 20:61098",
                1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File vcf = executeTest("testGATKVersionInVCF", spec).first.get(0);
        final VCFCodec codec = new VCFCodec();
        final VCFHeader header = (VCFHeader) codec.readActualHeader(codec.makeSourceFromStream(new FileInputStream(vcf)));

        // go through the metadata headers and look for ones that start with the GATK_COMMAND_LINE_KEY
        VCFHeaderLine versionLine = null;
        for ( final VCFHeaderLine headerLine : header.getMetaDataInInputOrder()) {
           if(headerLine.getKey().startsWith(GATKVCFUtils.GATK_COMMAND_LINE_KEY)) {
               versionLine = headerLine;
               break;
           }
        }
        Assert.assertNotNull(versionLine);
        Assert.assertTrue(versionLine.toString().contains("TestPrintVariantsWalker"));
    }

    @Test(enabled = true)
    public void testMultipleGATKVersionsInVCF() throws Exception {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintVariantsWalker -R " + b37KGReference +
                " -V " + privateTestDir + "gatkCommandLineInHeader.vcf"
                + " -o %s",
                1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File vcf = executeTest("testMultipleGATKVersionsInVCF", spec).first.get(0);
        final VCFCodec codec = new VCFCodec();
        final VCFHeader header = (VCFHeader) codec.readActualHeader(codec.makeSourceFromStream(new FileInputStream(vcf)));

        boolean foundHC = false;
        boolean foundPV = false;
        for ( final VCFHeaderLine line : header.getMetaDataInInputOrder() ) {
            if ( line.getKey().startsWith(GATKVCFUtils.GATK_COMMAND_LINE_KEY) ) {
                if ( line.toString().contains("HaplotypeCaller") ) {
                    Assert.assertFalse(foundHC);
                    foundHC = true;
                }
                if ( line.toString().contains("TestPrintVariantsWalker") ) {
                    Assert.assertFalse(foundPV);
                    foundPV = true;
                }
            }
        }

        Assert.assertTrue(foundHC, "Didn't find HaplotypeCaller command line header field");
        Assert.assertTrue(foundPV, "Didn't find TestPrintVariantsWalker command line header field");
    }

    @Test(enabled = true)
    public void testMultipleGATKVersionsSameWalkerInVCF() throws Exception {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintVariantsWalker -R " + b37KGReference +
                " -V " + privateTestDir + "gatkCommandLineExistsInHeader.vcf"
                + " -o %s",
                1, Arrays.asList(""));
        spec.disableShadowBCF();
        final File vcf = executeTest("testMultipleGATKVersionsSameWalkerInVCF", spec).first.get(0);
        final VCFCodec codec = new VCFCodec();
        final VCFHeader header = (VCFHeader) codec.readActualHeader(codec.makeSourceFromStream(new FileInputStream(vcf)));

        boolean foundFirstWalker  = false;
        boolean foundSecondWalker = false;
        for ( final VCFHeaderLine line : header.getMetaDataInInputOrder() ) {
            if ( line.getKey().startsWith(GATKVCFUtils.GATK_COMMAND_LINE_KEY) ) {
                // check if we found the second walker command line header field key
                if ( line.getKey().contains("TestPrintVariantsWalker.2") ) {
                    Assert.assertFalse(foundSecondWalker);
                    foundSecondWalker = true;
                }
                // otherwise if this is not the second walker command but contains the same
                // walker name, then it is the first occurrence.  If we somehow got more than
                // two occurrences of this walker, the Assert.assertFalse(foundFirstWalker);
                // will catch this
                else if ( line.getKey().contains("TestPrintVariantsWalker") ) {
                    Assert.assertFalse(foundFirstWalker);
                    foundFirstWalker = true;
                }
            }
        }

        Assert.assertTrue(foundFirstWalker, "Didn't find TestPrintVariantsWalker command line header field");
        Assert.assertTrue(foundSecondWalker, "Didn't find (second) TestPrintVariantsWalker command line header field");
    }

    // --------------------------------------------------------------------------------
    //
    // Test that defaultBaseQualities actually works
    //
    // --------------------------------------------------------------------------------

    public WalkerTestSpec testDefaultBaseQualities(final Integer value, final String md5) {
        return new WalkerTestSpec("-T TestPrintReadsWalker -R " + b37KGReference + " -I " + privateTestDir + "/baseQualitiesToFix.bam -o %s"
                + (value != null ? " --defaultBaseQualities " + value : ""),
                1, Arrays.asList(md5));
    }

    @Test()
    public void testDefaultBaseQualities20() {
        executeTest("testDefaultBaseQualities20", testDefaultBaseQualities(20, "90a450f74554bbd2cc3a9e0f9de68e26"));
    }

    @Test()
    public void testDefaultBaseQualities30() {
        executeTest("testDefaultBaseQualities30", testDefaultBaseQualities(30, "ec11db4173ce3b8e43997f00dab5ae26"));
    }

    @Test(expectedExceptions = Exception.class)
    public void testDefaultBaseQualitiesNoneProvided() {
        executeTest("testDefaultBaseQualitiesNoneProvided", testDefaultBaseQualities(null, ""));
    }

    // --------------------------------------------------------------------------------
    //
    // Test engine-level cigar consolidation
    //
    // --------------------------------------------------------------------------------

    @Test
    public void testGATKEngineConsolidatesCigars() {
        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintReadsWalker" +
                                                       " -R " + b37KGReference +
                                                       " -I " + privateTestDir + "zero_length_cigar_elements.bam" +
                                                       " -o %s",
                                                       1, Arrays.asList(""));  // No MD5s; we only want to check the cigar

        final File outputBam = executeTest("testGATKEngineConsolidatesCigars", spec).first.get(0);
        final SAMFileReader reader = new SAMFileReader(outputBam);
        reader.setValidationStringency(ValidationStringency.SILENT);

        final SAMRecord read = reader.iterator().next();
        reader.close();

        // Original cigar was 0M3M0M8M. Check that it's been consolidated after running through the GATK engine:
        Assert.assertEquals(read.getCigarString(), "11M", "Cigar 0M3M0M8M not consolidated correctly by the engine");
    }

    // --------------------------------------------------------------------------------
    //
    // Test on-the-fly sample renaming
    //
    // --------------------------------------------------------------------------------

    // On-the-fly sample renaming test case: one single-sample bam with multiple read groups
    @Test
    public void testOnTheFlySampleRenamingWithSingleBamFile() throws IOException {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam  myNewSampleName"));

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintReadsWalker" +
                                                       " -R " + b37KGReference +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam" +
                                                       " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                                                       " -o %s",
                                                       1, Arrays.asList(""));  // No MD5s; we only want to check the read groups

        final File outputBam = executeTest("testOnTheFlySampleRenamingWithSingleBamFile", spec).first.get(0);
        final SAMFileReader reader = new SAMFileReader(outputBam);

        for ( final SAMReadGroupRecord readGroup : reader.getFileHeader().getReadGroups() ) {
            Assert.assertEquals(readGroup.getSample(), "myNewSampleName", String.format("Sample for read group %s not renamed correctly", readGroup.getId()));
        }

        reader.close();
    }

    // On-the-fly sample renaming test case: three single-sample bams with multiple read groups per bam
    @Test
    public void testOnTheFlySampleRenamingWithMultipleBamFiles() throws IOException {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam  newSampleFor12878",
                              privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12891.HEADERONLY.bam  newSampleFor12891",
                              privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12892.HEADERONLY.bam  newSampleFor12892"));

        final Map<String, String> readGroupToNewSampleMap = new HashMap<>();
        for ( String inputBamID : Arrays.asList("12878", "12891", "12892") ) {
            final File inputBam = new File(privateTestDir + String.format("CEUTrio.HiSeq.WGS.b37.NA%s.HEADERONLY.bam", inputBamID));
            final SAMFileReader inputBamReader = new SAMFileReader(inputBam);
            final String newSampleName = String.format("newSampleFor%s", inputBamID);
            for ( final SAMReadGroupRecord readGroup : inputBamReader.getFileHeader().getReadGroups() ) {
                readGroupToNewSampleMap.put(readGroup.getId(), newSampleName);
            }
            inputBamReader.close();
        }

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintReadsWalker" +
                                                       " -R " + b37KGReference +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam" +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12891.HEADERONLY.bam" +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12892.HEADERONLY.bam" +
                                                       " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                                                       " -o %s",
                                                       1, Arrays.asList(""));  // No MD5s; we only want to check the read groups

        final File outputBam = executeTest("testOnTheFlySampleRenamingWithMultipleBamFiles", spec).first.get(0);
        final SAMFileReader outputBamReader = new SAMFileReader(outputBam);

        int totalReadGroupsSeen = 0;
        for ( final SAMReadGroupRecord readGroup : outputBamReader.getFileHeader().getReadGroups() ) {
            Assert.assertEquals(readGroup.getSample(), readGroupToNewSampleMap.get(readGroup.getId()),
                                String.format("Wrong sample for read group %s after on-the-fly renaming", readGroup.getId()));
            totalReadGroupsSeen++;
        }

        Assert.assertEquals(totalReadGroupsSeen, readGroupToNewSampleMap.size(), "Wrong number of read groups encountered in output bam file");

        outputBamReader.close();
    }

    // On-the-fly sample renaming test case: three single-sample bams with multiple read groups per bam,
    //                                       performing renaming in only SOME of the bams
    @Test
    public void testOnTheFlySampleRenamingWithMultipleBamFilesPartialRename() throws IOException {
        // Rename samples for NA12878 and NA12892, but not for NA12891
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam  newSampleFor12878",
                              privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12892.HEADERONLY.bam  newSampleFor12892"));

        final Map<String, String> readGroupToNewSampleMap = new HashMap<>();
        for ( String inputBamID : Arrays.asList("12878", "12891", "12892") ) {
            final File inputBam = new File(privateTestDir + String.format("CEUTrio.HiSeq.WGS.b37.NA%s.HEADERONLY.bam", inputBamID));
            final SAMFileReader inputBamReader = new SAMFileReader(inputBam);

            // Special-case NA12891, which we're not renaming:
            final String newSampleName = inputBamID.equals("12891") ? "NA12891" : String.format("newSampleFor%s", inputBamID);

            for ( final SAMReadGroupRecord readGroup : inputBamReader.getFileHeader().getReadGroups() ) {
                readGroupToNewSampleMap.put(readGroup.getId(), newSampleName);
            }
            inputBamReader.close();
        }

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintReadsWalker" +
                                                       " -R " + b37KGReference +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam" +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12891.HEADERONLY.bam" +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12892.HEADERONLY.bam" +
                                                       " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                                                       " -o %s",
                                                       1, Arrays.asList(""));  // No MD5s; we only want to check the read groups

        final File outputBam = executeTest("testOnTheFlySampleRenamingWithMultipleBamFilesPartialRename", spec).first.get(0);
        final SAMFileReader outputBamReader = new SAMFileReader(outputBam);

        int totalReadGroupsSeen = 0;
        for ( final SAMReadGroupRecord readGroup : outputBamReader.getFileHeader().getReadGroups() ) {
            Assert.assertEquals(readGroup.getSample(), readGroupToNewSampleMap.get(readGroup.getId()),
                                String.format("Wrong sample for read group %s after on-the-fly renaming", readGroup.getId()));
            totalReadGroupsSeen++;
        }

        Assert.assertEquals(totalReadGroupsSeen, readGroupToNewSampleMap.size(), "Wrong number of read groups encountered in output bam file");

        outputBamReader.close();
    }

    // On-the-fly sample renaming test case: two single-sample bams with read group collisions
    @Test
    public void testOnTheFlySampleRenamingWithReadGroupCollisions() throws IOException {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam  newSampleFor12878",
                              privateTestDir + "CEUTrio.HiSeq.WGS.b37.READ_GROUP_COLLISIONS_WITH_NA12878.HEADERONLY.bam  newSampleForNot12878"));

        final Set<String> na12878ReadGroups = new HashSet<>();
        final SAMFileReader inputBamReader = new SAMFileReader(new File(privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam"));
        for ( final SAMReadGroupRecord readGroup : inputBamReader.getFileHeader().getReadGroups() ) {
            na12878ReadGroups.add(readGroup.getId());
        }
        inputBamReader.close();

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintReadsWalker" +
                                                       " -R " + b37KGReference +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.HEADERONLY.bam" +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.READ_GROUP_COLLISIONS_WITH_NA12878.HEADERONLY.bam" +
                                                       " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                                                       " -o %s",
                                                       1, Arrays.asList(""));  // No MD5s; we only want to check the read groups

        final File outputBam = executeTest("testOnTheFlySampleRenamingWithReadGroupCollisions", spec).first.get(0);
        final SAMFileReader outputBamReader = new SAMFileReader(outputBam);

        int totalReadGroupsSeen = 0;
        for ( final SAMReadGroupRecord readGroup : outputBamReader.getFileHeader().getReadGroups() ) {
            String expectedSampleName = "";
            if ( na12878ReadGroups.contains(readGroup.getId()) ) {
                expectedSampleName = "newSampleFor12878";
            }
            else {
                expectedSampleName = "newSampleForNot12878";
            }

            Assert.assertEquals(readGroup.getSample(), expectedSampleName,
                                String.format("Wrong sample for read group %s after on-the-fly renaming", readGroup.getId()));
            totalReadGroupsSeen++;
        }

        Assert.assertEquals(totalReadGroupsSeen, na12878ReadGroups.size() * 2, "Wrong number of read groups encountered in output bam file");

        outputBamReader.close();
    }

    // On-the-fly sample renaming test case: a multi-sample bam (this should generate a UserException)
    @Test
    public void testOnTheFlySampleRenamingWithMultiSampleBam() throws IOException {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "CEUTrio.HiSeq.WGS.b37.MERGED.HEADERONLY.bam  myNewSampleName"));

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintReadsWalker" +
                                                       " -R " + b37KGReference +
                                                       " -I " + privateTestDir + "CEUTrio.HiSeq.WGS.b37.MERGED.HEADERONLY.bam" +
                                                       " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                                                       " -o %s",
                                                       1,
                                                       UserException.class); // expecting a UserException here

        executeTest("testOnTheFlySampleRenamingWithMultiSampleBam", spec);
    }

    // On-the-fly sample renaming test case: ensure that walkers can see the remapped sample names in individual reads
    @Test
    public void testOnTheFlySampleRenamingVerifyWalkerSeesNewSamplesInReads() throws IOException {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam  myNewSampleName"));

        final WalkerTestSpec spec = new WalkerTestSpec(" -T OnTheFlySampleRenamingVerifyingTestWalker" +
                                                       " -R " + b37KGReference +
                                                       " -I " + privateTestDir + "NA12878.HiSeq.b37.chr20.10_11mb.bam" +
                                                       " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                                                       " --newSampleName myNewSampleName" +
                                                       " -L 20:10000000-10001000",
                                                       1, Arrays.asList(""));

        // Test is a success if our custom walker doesn't throw an exception
        executeTest("testOnTheFlySampleRenamingVerifyWalkerSeesNewSamplesInReads", spec);
    }

    @Test
    public void testOnTheFlySampleRenamingSingleSampleVCF() throws IOException {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "NA12878.WGS.b37.chr20.firstMB.vcf  newSampleForNA12878"));

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintVariantsWalker" +
                " -R " + b37KGReference +
                " -V " + privateTestDir + "NA12878.WGS.b37.chr20.firstMB.vcf" +
                " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                " -o %s",
                1,
                Arrays.asList("")); // No MD5s -- we will inspect the output file manually

        final File outputVCF = executeTest("testOnTheFlySampleRenamingSingleSampleVCF", spec).first.get(0);
        verifySampleRenaming(outputVCF, "newSampleForNA12878");
    }

    private void verifySampleRenaming( final File outputVCF, final String newSampleName ) throws IOException {
        final Pair<VCFHeader, VCIterable<LineIterator>> headerAndVCIter = VCIterable.readAllVCs(outputVCF, new VCFCodec());
        final VCFHeader header = headerAndVCIter.getFirst();
        final VCIterable<LineIterator> iter = headerAndVCIter.getSecond();

        // Verify that sample renaming occurred at both the header and record levels (checking only the first 10 records):

        Assert.assertEquals(header.getGenotypeSamples().size(), 1, "Wrong number of samples in output vcf header");
        Assert.assertEquals(header.getGenotypeSamples().get(0), newSampleName, "Wrong sample name in output vcf header");

        int recordCount = 0;
        while ( iter.hasNext() && recordCount < 10 ) {
            final VariantContext vcfRecord = iter.next();
            Assert.assertEquals(vcfRecord.getSampleNames().size(), 1, "Wrong number of samples in output vcf record");
            Assert.assertEquals(vcfRecord.getSampleNames().iterator().next(), newSampleName, "Wrong sample name in output vcf record");
            recordCount++;
        }
    }

    @Test
    public void testOnTheFlySampleRenamingVerifyWalkerSeesNewSamplesInVCFRecords() throws Exception {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "samplerenametest_single_sample_gvcf.vcf    FOOSAMPLE"));

        final WalkerTestSpec spec = new WalkerTestSpec(" -T OnTheFlySampleRenamingVerifyingRodWalker" +
                " -R " + hg19Reference +
                " -V " + privateTestDir + "samplerenametest_single_sample_gvcf.vcf" +
                " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                " --expectedSampleName FOOSAMPLE" +
                " -o %s",
                1,
                Arrays.asList("")); // No MD5s -- custom walker will throw an exception if there's a problem

        executeTest("testOnTheFlySampleRenamingVerifyWalkerSeesNewSamplesInVCFRecords", spec);
    }

    @Test
    public void testOnTheFlySampleRenamingMultiSampleVCF() throws Exception {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "vcf/vcfWithGenotypes.vcf  badSample"));

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintVariantsWalker" +
                " -R " + b37KGReference +
                " -V " + privateTestDir + "vcf/vcfWithGenotypes.vcf" +
                " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                " -o %s",
                1,
                UserException.class); // expecting a UserException here

        executeTest("testOnTheFlySampleRenamingMultiSampleVCF", spec);
    }

    @Test
    public void testOnTheFlySampleRenamingSitesOnlyVCF() throws Exception {
        final File sampleRenameMapFile = createTestSampleRenameMapFile(
                Arrays.asList(privateTestDir + "vcf/vcfWithoutGenotypes.vcf  badSample"));

        final WalkerTestSpec spec = new WalkerTestSpec(" -T TestPrintVariantsWalker" +
                " -R " + b37KGReference +
                " -V " + privateTestDir + "vcf/vcfWithoutGenotypes.vcf" +
                " --sample_rename_mapping_file " + sampleRenameMapFile.getAbsolutePath() +
                " -o %s",
                1,
                UserException.class); // expecting a UserException here

        executeTest("testOnTheFlySampleRenamingSitesOnlyVCF", spec);
    }

    private File createTestSampleRenameMapFile( final List<String> contents ) throws IOException {
        final File mapFile = createTempFile("TestSampleRenameMapFile", ".tmp");
        final PrintWriter writer = new PrintWriter(mapFile);

        for ( final String line : contents ) {
            writer.println(line);
        }
        writer.close();

        return mapFile;
    }

    public static class OnTheFlySampleRenamingVerifyingTestWalker extends ReadWalker<Integer, Integer> {
        @Argument(fullName = "newSampleName", shortName = "newSampleName", doc = "", required = true)
        String newSampleName = null;

        public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
            if ( ! newSampleName.equals(read.getReadGroup().getSample()) ) {
                throw new IllegalStateException(String.format("Encountered read with the wrong sample name. Expected %s found %s",
                                                              newSampleName, read.getReadGroup().getSample()));
            }

            return 1;
        }

        public Integer reduceInit() { return 0; }
        public Integer reduce(Integer value, Integer sum) { return value + sum; }
    }

    public static class OnTheFlySampleRenamingVerifyingRodWalker extends RodWalker<Integer, Integer> {
        @Argument(fullName = "expectedSampleName", shortName = "expectedSampleName", doc = "", required = true)
        String expectedSampleName = null;

        @Output
        PrintStream out;

        @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
        public RodBinding<VariantContext> variants;

        public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
            if ( tracker == null ) {
                return 0;
            }

            for ( final VariantContext vc : tracker.getValues(variants, context.getLocation()) ) {
                if ( vc.getSampleNames().size() != 1 ) {
                    throw new IllegalStateException("Encountered a vcf record with num samples != 1");
                }

                final String actualSampleName = vc.getSampleNames().iterator().next();
                if ( ! expectedSampleName.equals(actualSampleName)) {
                    throw new IllegalStateException(String.format("Encountered vcf record with wrong sample name. Expected %s found %s",
                                                                  expectedSampleName, actualSampleName));
                }
            }

            return 1;
        }

        public Integer reduceInit() {
            return 0;
        }

        public Integer reduce(Integer counter, Integer sum) {
            return counter + sum;
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Test output file-specific options
    //
    // --------------------------------------------------------------------------------

    //Returns the output file
    private File testBAMFeatures(final String args, final String md5) {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintReadsWalker -R " + b37KGReference +
                " -I " + privateTestDir + "NA20313.highCoverageRegion.bam"
                + " --no_pg_tag -o %s " + args,
                1, Arrays.asList(".bam"), Arrays.asList(md5));
        return executeTest("testBAMFeatures: "+args, spec).first.get(0);
    }

    @Test
    public void testSAMWriterFeatures() {
        testBAMFeatures("-compress 0", "49228d4f5b14c4cfed4a09372eb71139");
        testBAMFeatures("-compress 9", "bc61a1b2b53a2ec7c63b533fa2f8701b");
        testBAMFeatures("-simplifyBAM", "a1127bab46674b165496b79bb9fa7964");

        //Validate MD5
        final String expectedMD5 = "c58b9114fc15b53655f2c03c819c29fd";
        final File md5Target = testBAMFeatures("--generate_md5", expectedMD5);
        final File md5File = new File(md5Target.getAbsoluteFile() + ".md5");
        md5File.deleteOnExit();
        Assert.assertTrue(md5File.exists(), "MD5 wasn't created");
        try {
            String md5 = new BufferedReader(new FileReader(md5File)).readLine();
            Assert.assertEquals(md5, expectedMD5, "Generated MD5 doesn't match expected");
        } catch (IOException e) {
            Assert.fail("Can't parse MD5 file", e);
        }

        //Validate that index isn't created
        final String unindexedBAM = testBAMFeatures("--disable_bam_indexing", expectedMD5).getAbsolutePath();
        Assert.assertTrue(!(new File(unindexedBAM+".bai").exists()) &&
                          !(new File(unindexedBAM.replace(".bam", ".bai")).exists()),
                          "BAM index was created even though it was disabled");
    }

    @DataProvider(name = "vcfFeaturesData")
    public Object[][] getVCFFeaturesData() {
        return new Object[][]{
                {"--sites_only", "6ef742ee6d9bcbc7b23f928c0e8a1d0e"},
                {"--bcf", "285549ca1a719a09fa95cfa129520621"}
        };
    }

    @Test(dataProvider = "vcfFeaturesData")
    public void testVCFFeatures(final String args, final String md5) {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintVariantsWalker -R " + b37KGReference +
                " -V " + privateTestDir + "CEUtrioTest.vcf"
                + " --no_cmdline_in_header -o %s " + args,
                1, Arrays.asList(md5));
        executeTest("testVCFFeatures: "+args, spec);
    }

    @DataProvider(name = "vcfFormatHandlingData")
    public Object[][] getVCFFormatHandlingData() {
        return new Object[][]{
                {true, "870f39e19ec89c8a09f7eca0f5c4bcb9"},
                {false, "baf9a1755d3b4e0ed25b03233e99ca91"}
        };
    }

    @Test(dataProvider = "vcfFormatHandlingData")
    public void testVCFFormatHandling(final boolean writeFullFormat, final String md5) {
        WalkerTestSpec spec = new WalkerTestSpec("-T TestPrintVariantsWalker -R " + b37KGReference +
                " -V " + privateTestDir + "ILLUMINA.wex.broad_phase2_baseline.20111114.both.exome.genotypes.1000.vcf"
                + " --no_cmdline_in_header -o %s "
                + " --fullyDecode " //Without this parameter, the FORMAT fields will be emitted unchanged.  Oops
                + (writeFullFormat ? "-writeFullFormat" : "") ,
                1, Arrays.asList(md5));
        executeTest("testVCFFormatHandling: "+(writeFullFormat ? "Untrimmed" : "Trimmed"), spec);
    }
}
