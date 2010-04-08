package org.broadinstitute.sting.gatk.datasources.shards;

import org.junit.Test;
import org.junit.Before;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
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
 *         Class ReadIntervalShardStrategyUnitTest
 *         <p/>
 *         Tests the ReadIntervalShardStrategy class
 *
 *
 * TODO: this test has been changed, since we now force interval shards to not subdivide if they're large,
 * so you should always get one shard back, reguardless of the specified length.  If this behavior changes
 * back in the future, the expected number of shards should change from 1 to > 1.
 *
 */
public class IntervalShardStrategyUnitTest extends BaseTest {

    private GenomeLocSortedSet mSortedSet = null;
    private SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    @Before
    public void setup() {
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
        mSortedSet = new GenomeLocSortedSet();
    }

    @Test
    public void testNoExceptionOnEmpty() {
        IntervalShardStrategy strat = new IntervalShardStrategy(100, mSortedSet,Shard.ShardType.LOCUS_INTERVAL);
    }

    @Test
    public void testSingleChromosomeFunctionality() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, 1, 1000);
        mSortedSet.add(loc);
        IntervalShardStrategy strat = new IntervalShardStrategy(100, mSortedSet,Shard.ShardType.LOCUS_INTERVAL);
        int counter = 0;
        Shard d = null;
        while (strat.hasNext()) {
            d = strat.next();
            counter++;
        }
        assertTrue(d instanceof IntervalShard);
        assertEquals(1, counter);
    }

    @Test
    public void testMultipleChromosomeFunctionality() {
        for (int x = 0; x < 5; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(x, 1, 1000);
            mSortedSet.add(loc);
        }
        IntervalShardStrategy strat = new IntervalShardStrategy(100, mSortedSet,Shard.ShardType.LOCUS_INTERVAL);
        int counter = 0;
        Shard d = null;
        while (strat.hasNext()) {
            d = strat.next();
            counter++;
        }
        assertTrue(d instanceof IntervalShard);
        assertEquals(5, counter);
    }

    @Test
    public void testOddSizeShardFunctionality() {
        for (int x = 0; x < 5; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(x, 1, 1000);
            mSortedSet.add(loc);
        }
        IntervalShardStrategy strat = new IntervalShardStrategy(789, mSortedSet,Shard.ShardType.LOCUS_INTERVAL);
        int counter = 0;
        while (strat.hasNext()) {
            Shard d = strat.next();
            assertEquals(1,d.getGenomeLocs().size());
            assertEquals(1, d.getGenomeLocs().get(0).getStart());
            assertEquals(1000, d.getGenomeLocs().get(0).getStop());
            counter++;
        }
        assertEquals(5, counter);
    }


    @Test
    public void testInfiniteShardSize() {
        for (int x = 0; x < 5; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(x, 1, 1000);
            mSortedSet.add(loc);
        }
        IntervalShardStrategy strat = new IntervalShardStrategy(Long.MAX_VALUE, mSortedSet,Shard.ShardType.LOCUS_INTERVAL);
        int counter = 0;
        while (strat.hasNext()) {
            Shard d = strat.next();
            assertEquals(1,d.getGenomeLocs().size());
            assertEquals(1000, d.getGenomeLocs().get(0).getStop());
            counter++;
        }
        assertEquals(5, counter);
    }

    @Test(expected = UnsupportedOperationException.class)
    public void testRemove() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, 1, 1000);
        mSortedSet.add(loc);
        IntervalShardStrategy strat = new IntervalShardStrategy(100, mSortedSet,Shard.ShardType.LOCUS_INTERVAL);
        strat.remove();
    }

}
