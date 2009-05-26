package org.broadinstitute.sting.gatk.dataSources.shards;

import org.junit.Test;
import org.junit.Before;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.sam.ArtificialSamUtils;
import org.broadinstitute.sting.BaseTest;
import net.sf.samtools.SAMFileHeader;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class ReadIntervalShardStrategyTest
 *         <p/>
 *         Tests the ReadIntervalShardStrategy class
 */
public class ReadIntervalShardStrategyTest extends BaseTest {

    private GenomeLocSortedSet mSortedSet = null;
    private SAMFileHeader header = ArtificialSamUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    @Before
    public void setup() {
        GenomeLoc.setupRefContigOrdering(header.getSequenceDictionary());
        mSortedSet = new GenomeLocSortedSet();
    }

    @Test(expected = StingException.class)
    public void testExceptionOnEmpty() {
        ReadIntervalShardStrategy strat = new ReadIntervalShardStrategy(header.getSequenceDictionary(), 100, mSortedSet);
    }

    @Test
    public void testSingleChromosomeFunctionality() {
        GenomeLoc loc = new GenomeLoc(1, 1, 1000);
        mSortedSet.add(loc);
        ReadIntervalShardStrategy strat = new ReadIntervalShardStrategy(header.getSequenceDictionary(), 100, mSortedSet);
        int counter = 0;
        while (strat.hasNext()) {
            Shard d = strat.next();
            counter++;
        }
        assertEquals(10, counter);
    }

    @Test
    public void testMultipleChromosomeFunctionality() {
        for (int x = 0; x < 5; x++) {
            GenomeLoc loc = new GenomeLoc(x, 1, 1000);
            mSortedSet.add(loc);
        }
        ReadIntervalShardStrategy strat = new ReadIntervalShardStrategy(header.getSequenceDictionary(), 100, mSortedSet);
        int counter = 0;
        while (strat.hasNext()) {
            Shard d = strat.next();
            counter++;
        }
        assertEquals(50, counter);
    }

    @Test
    public void testOddSizeShardFunctionality() {
        for (int x = 0; x < 5; x++) {
            GenomeLoc loc = new GenomeLoc(x, 1, 1000);
            mSortedSet.add(loc);
        }
        ReadIntervalShardStrategy strat = new ReadIntervalShardStrategy(header.getSequenceDictionary(), 789, mSortedSet);
        int counter = 0;
        while (strat.hasNext()) {
            Shard d = strat.next();
            if (counter % 2 == 0) {
                assertEquals(1, d.getGenomeLoc().getStart());
                assertEquals(789, d.getGenomeLoc().getStop());
            } else {
                assertEquals(790, d.getGenomeLoc().getStart());
                assertEquals(1000, d.getGenomeLoc().getStop());
            }
            counter++;
        }
        assertEquals(10, counter);
    }

    @Test(expected = UnsupportedOperationException.class)
    public void testRemove() {
        GenomeLoc loc = new GenomeLoc(1, 1, 1000);
        mSortedSet.add(loc);
        ReadIntervalShardStrategy strat = new ReadIntervalShardStrategy(header.getSequenceDictionary(), 100, mSortedSet);
        strat.remove();
    }

}
