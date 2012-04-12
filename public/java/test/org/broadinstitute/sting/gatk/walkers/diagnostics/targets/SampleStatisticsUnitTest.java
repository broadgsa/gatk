/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SampleStatisticsUnitTest/* extends BaseTest */ {

    @DataProvider(name = "QuartileValues")
    public Object[][] getQuantileValues() {

        int[] a1 = {5};
        int[] a2 = {1, 2};
        int[] a5 = {10, 20, 30, 40, 50};
        int[] a10 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};


        return new Object[][]{
                new Object[]{a1, 0.5, 5},
                new Object[]{a1, 0, 5},
                new Object[]{a1, 1, 5},
                new Object[]{a2, 0.5, 1.5},
                new Object[]{a2, 0.25, 1},
                new Object[]{a2, 0.75, 2},
                new Object[]{a5, 0.5, 30},
                new Object[]{a5, 0.25, 20},
                new Object[]{a5, 0.75, 40},
                new Object[]{a5, 0, -1},
                new Object[]{a10, 0.5, 5.5},
                new Object[]{a10, 0.25, 3},
                new Object[]{a10, 0.75, 8}
        };
    }

    @Test(dataProvider = "QuartileValues")
    public void testGetQuartile(int[] dataList, double percentage, double expected) {
        Assert.assertEquals(SampleStatistics.getQuartile(dataList, percentage), expected);

    }

    @DataProvider(name = "ReadsAndMates")
    public Object[][] getReadAndMates() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

        GATKSAMRecord noPair = ArtificialSAMUtils.createArtificialRead(header, "test", 0, 100, 50);
        GATKSAMRecord good = ArtificialSAMUtils.createPair(header, "test", 30, 100, 150, true, false).get(0);
        GATKSAMRecord bigInsertSize = ArtificialSAMUtils.createPair(header, "test", 30, 100, 151, true, false).get(0);
        GATKSAMRecord inverted = ArtificialSAMUtils.createPair(header, "test", 30, 151, 150, true, false).get(0);
        GATKSAMRecord sameOrientation = ArtificialSAMUtils.createPair(header, "test", 30, 100, 151, true, true).get(0);

        GATKSAMRecord pairNotMapped = ArtificialSAMUtils.createPair(header, "test", 30, 100, 140, true, false).get(1);
        pairNotMapped.setMateUnmappedFlag(true);

        // finish test
        return new Object[][]{
                new Object[]{noPair, false},
                new Object[]{good, true},
                new Object[]{bigInsertSize, false},
                new Object[]{inverted, false},
                new Object[]{sameOrientation, false},
                new Object[]{pairNotMapped, false}
        };
    }

    @Test(dataProvider = "ReadsAndMates")
    public void testHasValidMate(GATKSAMRecord read, boolean expected) {
        //50 is out maximum insert size
        Assert.assertEquals(new SampleStatistics(GenomeLoc.UNMAPPED).hasValidMate(read, ThresHolder.DEFAULTS), expected);
    }

}
