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

package org.broadinstitute.gatk.engine.datasources.rmd;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.refdata.utils.*;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class ReferenceOrderedQueryDataPoolUnitTest extends BaseTest{
    @Test
    public void testCloseFilePointers() throws IOException {
        // Build up query parameters
        File file = new File(BaseTest.privateTestDir + "NA12878.hg19.example1.vcf");
        RMDTriplet triplet = new RMDTriplet("test", "VCF", file.getAbsolutePath(), RMDTriplet.RMDStorageType.FILE, new Tags());
        IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(BaseTest.hg19Reference));
        GenomeLocParser parser = new GenomeLocParser(seq);
        GenomeLoc loc = parser.createGenomeLoc("20", 1, 100000);
        TestRMDTrackBuilder builder = new TestRMDTrackBuilder(seq.getSequenceDictionary(), parser);

        // Create the query data pool
        ReferenceOrderedQueryDataPool pool = new ReferenceOrderedQueryDataPool(triplet, builder, seq.getSequenceDictionary(), parser);

        for (int i = 0; i < 3; i++) {
            // Ensure our tribble iterators are closed.
            CheckableCloseableTribbleIterator.clearThreadIterators();
            Assert.assertTrue(CheckableCloseableTribbleIterator.getThreadIterators().isEmpty(), "Tribble iterators list was not cleared.");

            // Request the the rodIterator
            LocationAwareSeekableRODIterator rodIterator = pool.iterator(new MappedStreamSegment(loc));

            // Run normal iteration over rodIterator
            Assert.assertTrue(rodIterator.hasNext(), "Rod iterator does not have a next value.");
            GenomeLoc rodIteratorLocation = rodIterator.next().getLocation();
            Assert.assertEquals(rodIteratorLocation.getContig(), "20", "Instead of chr 20 rod iterator was at location " + rodIteratorLocation);

            // Check that the underlying tribbleIterators are still open.
            List<CheckableCloseableTribbleIterator<? extends Feature>> tribbleIterators = CheckableCloseableTribbleIterator.getThreadIterators();
            Assert.assertFalse(tribbleIterators.isEmpty(), "Tribble iterators list is empty");
            for (CheckableCloseableTribbleIterator<? extends Feature> tribbleIterator: tribbleIterators) {
                Assert.assertFalse(tribbleIterator.isClosed(), "Tribble iterator is closed but should be still open.");
            }

            // Releasing the rodIterator should close the underlying tribbleIterator.
            pool.release(rodIterator);

            // Check that the underlying tribbleIterators are now closed.
            for (CheckableCloseableTribbleIterator<? extends Feature> tribbleIterator: tribbleIterators) {
                Assert.assertTrue(tribbleIterator.isClosed(), "Tribble iterator is open but should be now closed.");
            }
        }

        // Extra cleanup.
        CheckableCloseableTribbleIterator.clearThreadIterators();
    }
}
