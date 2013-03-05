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

package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.commandline.ArgumentException;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.readutils.PrintReads;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Tests selected functionality in the GenomeAnalysisEngine class
 */
public class GenomeAnalysisEngineUnitTest extends BaseTest {

    @Test(expectedExceptions=ArgumentException.class)
    public void testDuplicateSamFileHandlingSingleDuplicate() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"), new Tags()));

        testEngine.setSAMFileIDs(samFiles);
        testEngine.checkForDuplicateSamFiles();
    }

    @Test(expectedExceptions=ArgumentException.class)
    public void testDuplicateSamFileHandlingMultipleDuplicates() throws Exception {
        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleNORG.bam"), new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleBAM.bam"),  new Tags()));
        samFiles.add(new SAMReaderID(new File("public/testdata/exampleNORG.bam"), new Tags()));

        testEngine.setSAMFileIDs(samFiles);
        testEngine.checkForDuplicateSamFiles();
    }

    @Test
    public void testEmptyIntervalSetHandling() throws Exception {
        GenomeLocParser genomeLocParser = new GenomeLocParser(ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000).getSequenceDictionary());

        GenomeAnalysisEngine testEngine = new GenomeAnalysisEngine();

        testEngine.setWalker(new PrintReads());
        testEngine.setIntervals(new GenomeLocSortedSet(genomeLocParser));

        testEngine.validateSuppliedIntervals();
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
