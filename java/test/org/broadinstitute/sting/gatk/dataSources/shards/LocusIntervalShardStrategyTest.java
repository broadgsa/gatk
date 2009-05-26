package org.broadinstitute.sting.gatk.dataSources.shards;

import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.ArtificialSamUtils;
import org.broadinstitute.sting.BaseTest;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertTrue;
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
 *         Class LocusIntervalShardStrategyTest
 *         <p/>
 *         Tests the LocusIntervalShardStrategy class.
 */
public class LocusIntervalShardStrategyTest extends BaseTest {
    private GenomeLocSortedSet mSortedSet = null;
    private SAMFileHeader header = ArtificialSamUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;
    private LocusIntervalShardStrategy strat = null;

    @Before
    public void setup() {
        GenomeLoc.setupRefContigOrdering(header.getSequenceDictionary());
        mSortedSet = new GenomeLocSortedSet();
    }

    @Test
    public void testOneToOneness() {
        for (int x = 0; x < 100; x++) {
            GenomeLoc loc = new GenomeLoc(0,(x*10)+1, (x*10)+8);
            mSortedSet.add(loc);
        }
        strat = new LocusIntervalShardStrategy(header.getSequenceDictionary(),mSortedSet);
        int counter = 0;
        while (strat.hasNext()) {
            ++counter;
            GenomeLoc loc = strat.next().getGenomeLoc();
            long stop = loc.getStop();
            long start = loc.getStart();
            long length =  stop - start;
            assertTrue(length == 7);
        }
        assertTrue(counter == 100);

    }

}
