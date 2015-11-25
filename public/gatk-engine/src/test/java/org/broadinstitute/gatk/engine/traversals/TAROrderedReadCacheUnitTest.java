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

package org.broadinstitute.gatk.engine.traversals;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.ArtificialBAMBuilder;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TAROrderedReadCacheUnitTest extends BaseTest {
    // example fasta index file, can be deleted if you don't use the reference
    private IndexedFastaSequenceFile seq;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
    }

    @DataProvider(name = "ReadCacheTestData")
    public Object[][] makeReadCacheTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int nReadsPerLocus : Arrays.asList(0, 1, 10, 100) ) {
            for ( final int nLoci : Arrays.asList(1, 10, 100) ) {
                for ( final int max : Arrays.asList(10, 50, 1000) ) {
                    for ( final boolean addAllAtOnce : Arrays.asList(true, false) ) {
                        tests.add(new Object[]{nReadsPerLocus, nLoci, max, addAllAtOnce});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadCacheTestData")
    public void testReadCache(final int nReadsPerLocus, final int nLoci, final int max, final boolean addAllAtOnce) {
        final TAROrderedReadCache cache = new TAROrderedReadCache(max);

        Assert.assertEquals(cache.getMaxCapacity(), max);
        Assert.assertEquals(cache.getNumDiscarded(), 0);
        Assert.assertEquals(cache.size(), 0);

        final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(seq, nReadsPerLocus, nLoci);
        final List<GATKSAMRecord> reads = bamBuilder.makeReads();

        if ( addAllAtOnce ) {
            cache.addAll(reads);
        } else {
            for ( final GATKSAMRecord read : reads ) {
                cache.add(read);
            }
        }

        final int nTotalReads = reads.size();
        final int nExpectedToKeep = Math.min(nTotalReads, max);
        final int nExpectedToDiscard = nTotalReads - nExpectedToKeep;
        Assert.assertEquals(cache.getNumDiscarded(), nExpectedToDiscard, "wrong number of reads discarded");
        Assert.assertEquals(cache.size(), nExpectedToKeep, "wrong number of reads kept");

        final List<GATKSAMRecord> cacheReads = cache.popCurrentReads();
        Assert.assertEquals(cache.size(), 0, "Should be no reads left");
        Assert.assertEquals(cache.getNumDiscarded(), 0, "should have reset stats");
        Assert.assertEquals(cacheReads.size(), nExpectedToKeep, "should have 1 read for every read we expected to keep");

        verifySortednessOfReads(cacheReads);
    }

    private void verifySortednessOfReads( final List<GATKSAMRecord> reads) {
        int lastStart = -1;
        for ( GATKSAMRecord read : reads ) {
            Assert.assertTrue(lastStart <= read.getAlignmentStart(), "Reads should be sorted but weren't.  Found read with start " + read.getAlignmentStart() + " while last was " + lastStart);
            lastStart = read.getAlignmentStart();
        }
    }
}
