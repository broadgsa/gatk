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

import org.broadinstitute.gatk.engine.walkers.TestCountReadsWalker;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Tests selected functionality in the GenomeAnalysisEngine class
 */
public class GenomeAnalysisEngineUnitTest extends BaseTest {

    @Test(expectedExceptions=UserException.class)
    public void testEmptySamFileListHandling() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();
        testEngine.setWalker(new TestCountReadsWalker()); //generalizable to any walker requiring reads

        //supply command line args so validateSuppliedReads() knows whether reads were passed in
        GATKArgumentCollection testArgs = new GATKArgumentCollection();
        testArgs.samFiles.add("empty.list");
        testEngine.setArguments(testArgs);

        //represents the empty list of samFiles read in from empty.list by CommandLineExecutable
        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();

        testEngine.setSAMFileIDs(samFiles);
        testEngine.validateSuppliedReads();
    }

    @Test(expectedExceptions=UserException.class)
    public void testDuplicateSamFileHandlingSingleDuplicate() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        samFiles.add(new SAMReaderID(new File(publicTestDir + "exampleBAM.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File(publicTestDir + "exampleBAM.bam"), new Tags()));

        testEngine.setSAMFileIDs(samFiles);
        testEngine.checkForDuplicateSamFiles();
    }

    @Test(expectedExceptions=UserException.class)
    public void testDuplicateSamFileHandlingMultipleDuplicates() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        samFiles.add(new SAMReaderID(new File(publicTestDir + "exampleBAM.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File(publicTestDir + "exampleNORG.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File(publicTestDir + "exampleBAM.bam"),  new Tags()));
        samFiles.add(new SAMReaderID(new File(publicTestDir + "exampleNORG.bam"), new Tags()));

        testEngine.setSAMFileIDs(samFiles);
        testEngine.checkForDuplicateSamFiles();
    }

    @Test(expectedExceptions=UserException.class)
    public void testDuplicateSamFileHandlingAbsoluteVsRelativePath() {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        final File relativePathToBAMFile = new File(publicTestDir + "exampleBAM.bam");
        final File absolutePathToBAMFile = new File(relativePathToBAMFile.getAbsolutePath());
        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        samFiles.add(new SAMReaderID(relativePathToBAMFile, new Tags()));
        samFiles.add(new SAMReaderID(absolutePathToBAMFile, new Tags()));

        testEngine.setSAMFileIDs(samFiles);
        testEngine.checkForDuplicateSamFiles();
    }

    @Test
    public void testEmptyIntervalSetHandling() throws Exception {
        GenomeLocParser genomeLocParser = new GenomeLocParser(ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000).getSequenceDictionary());

        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        testEngine.setWalker(new TestCountReadsWalker());
        testEngine.setIntervals(new GenomeLocSortedSet(genomeLocParser));

        testEngine.validateSuppliedIntervals();
    }

    @Test
    public void testLoadWellFormedSampleRenameMapFile() throws IOException {
        final File mapFile = createTestSampleRenameMapFile(Arrays.asList("/foo/bar/first.bam    newSample1",
                                                                         "/foo/bar/second.bam        newSample2",
                                                                         "/foo/bar2/third.bam newSample3",
                                                                         "/foo/bar2/fourth.bam new sample    4",
                                                                         "/foo/bar2/fifth.bam     new   sample     5    "));
        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        final Map<String, String> renameMap = engine.loadSampleRenameMap(mapFile);

        Assert.assertEquals(renameMap.size(), 5, "Sample rename map was wrong size after loading from file");

        final Iterator<String> expectedResultsIterator = Arrays.asList(
                        "/foo/bar/first.bam",   "newSample1", 
                        "/foo/bar/second.bam",  "newSample2", 
                        "/foo/bar2/third.bam",  "newSample3",
                        "/foo/bar2/fourth.bam", "new sample    4",
                        "/foo/bar2/fifth.bam",  "new   sample     5"
        ).iterator();
        while ( expectedResultsIterator.hasNext() ) {
            final String expectedKey = expectedResultsIterator.next();
            final String expectedValue = expectedResultsIterator.next();

            Assert.assertNotNull(renameMap.get(expectedKey), String.format("Entry for %s not found in sample rename map", expectedKey));
            Assert.assertEquals(renameMap.get(expectedKey), expectedValue, "Wrong value in sample rename map for " + expectedKey);
        }
    }

    @DataProvider(name = "MalformedSampleRenameMapFileDataProvider")
    public Object[][] generateMalformedSampleRenameMapFiles() throws IOException {
        final List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{"testLoadSampleRenameMapFileNonExistentFile",
                               new File("/foo/bar/nonexistent")});
        tests.add(new Object[]{"testLoadSampleRenameMapFileMalformedLine",
                               createTestSampleRenameMapFile(Arrays.asList("/path/to/foo.bam"))});
        tests.add(new Object[]{"testLoadSampleRenameMapFileNonAbsoluteBamPath",
                               createTestSampleRenameMapFile(Arrays.asList("relative/path/to/foo.bam newSample"))});
        tests.add(new Object[]{"testLoadSampleRenameMapFileDuplicateBamPath",
                               createTestSampleRenameMapFile(Arrays.asList("/path/to/dupe.bam newSample1",
                                                                           "/path/to/dupe.bam newSample2"))});
        tests.add(new Object[]{"testLoadSampleRenameMapFileTabInSampleName",
                               createTestSampleRenameMapFile(Arrays.asList("/path/to/stuff.bam some wonky\tsample   "))});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MalformedSampleRenameMapFileDataProvider", expectedExceptions = UserException.class)
    public void testLoadMalformedSampleRenameMapFile( final String testName, final File mapFile ) {
        logger.info("Executing test " + testName);

        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        final Map<String, String> renameMap = engine.loadSampleRenameMap(mapFile);
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

    ///////////////////////////////////////////////////
    // Test the ReadTransformer ordering enforcement //
    ///////////////////////////////////////////////////

    public static class TestReadTransformer extends ReadTransformer {

        private OrderingConstraint orderingConstraint = OrderingConstraint.DO_NOT_CARE;
        private boolean enabled;

        protected TestReadTransformer(final OrderingConstraint orderingConstraint) {
            this.orderingConstraint = orderingConstraint;
            enabled = true;
        }

        // need this because PackageUtils will pick up this class as a possible ReadTransformer
        protected TestReadTransformer() {
            enabled = false;
        }

        @Override
        public OrderingConstraint getOrderingConstraint() { return orderingConstraint; }

        @Override
        public ApplicationTime initializeSub(final GenomeAnalysisEngine engine, final Walker walker) { return ApplicationTime.HANDLED_IN_WALKER; }

        @Override
        public boolean enabled() { return enabled; }

        @Override
        public GATKSAMRecord apply(final GATKSAMRecord read) { return read; }

    }

    @DataProvider(name = "ReadTransformerData")
    public Object[][] makeReadTransformerData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final ReadTransformer.OrderingConstraint orderingConstraint1 : ReadTransformer.OrderingConstraint.values() ) {
            for ( final ReadTransformer.OrderingConstraint orderingConstraint2 : ReadTransformer.OrderingConstraint.values() ) {
                for ( final ReadTransformer.OrderingConstraint orderingConstraint3 : ReadTransformer.OrderingConstraint.values() ) {
                    tests.add(new Object[]{orderingConstraint1, orderingConstraint2, orderingConstraint3});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadTransformerData")
    public void testReadTransformer(final ReadTransformer.OrderingConstraint oc1, final ReadTransformer.OrderingConstraint oc2, final ReadTransformer.OrderingConstraint oc3) {

        final GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();
        final List<ReadTransformer> readTransformers = new ArrayList<ReadTransformer>(3);
        readTransformers.add(new TestReadTransformer(oc1));
        readTransformers.add(new TestReadTransformer(oc2));
        readTransformers.add(new TestReadTransformer(oc3));

        final boolean shouldThrowException = numWithConstraint(ReadTransformer.OrderingConstraint.MUST_BE_FIRST, oc1, oc2, oc3) > 1 ||
                numWithConstraint(ReadTransformer.OrderingConstraint.MUST_BE_LAST, oc1, oc2, oc3) > 1;

        try {
            testEngine.setReadTransformers(readTransformers);

            Assert.assertFalse(shouldThrowException);
            Assert.assertEquals(testEngine.getReadTransformers().size(), 3);

            Assert.assertTrue(testEngine.getReadTransformers().get(1).getOrderingConstraint() != ReadTransformer.OrderingConstraint.MUST_BE_FIRST);
            Assert.assertTrue(testEngine.getReadTransformers().get(2).getOrderingConstraint() != ReadTransformer.OrderingConstraint.MUST_BE_FIRST);
            Assert.assertTrue(testEngine.getReadTransformers().get(0).getOrderingConstraint() != ReadTransformer.OrderingConstraint.MUST_BE_LAST);
            Assert.assertTrue(testEngine.getReadTransformers().get(1).getOrderingConstraint() != ReadTransformer.OrderingConstraint.MUST_BE_LAST);
        } catch (UserException.IncompatibleReadFiltersException e) {
            Assert.assertTrue(shouldThrowException);
        }
    }

    private int numWithConstraint(final ReadTransformer.OrderingConstraint target, final ReadTransformer.OrderingConstraint... constraints ) {
        int count = 0;
        for ( final ReadTransformer.OrderingConstraint constraint : constraints ) {
            if ( constraint == target )
                count++;
        }
        return count;
    }
}
