package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategyFactory;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;


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
 *         Class SAMByIntervalUnitTest
 *         <p/>
 *         Test that the SAM data source behaves well given intervals
 */
public class SAMByIntervalUnitTest extends BaseTest {
    private List<File> fl;
    ShardStrategy shardStrategy;
    Reads reads;
    private int targetReadCount = 14;


    // constants we use throughout the tests
    protected final int READ_COUNT;
    protected final int ENDING_CHROMO;
    protected final int STARTING_CHROMO;
    protected final int UNMAPPED_READ_COUNT;

    public SAMByIntervalUnitTest() {
        READ_COUNT = 100;
        ENDING_CHROMO = 10;
        STARTING_CHROMO = 1;
        UNMAPPED_READ_COUNT = 1000;
    }

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {

        fl = new ArrayList<File>();

        // sequence
        //seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly17/v0/Homo_sapiens_assembly17.fasta"));
        //GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());

        // setup the test files
        fl.add(new File(validationDataLocation + "index_test.bam"));
        reads = new Reads(fl);
    }


    /** run a test on data over a specific interval */
    private void testRead( int start, int stop, int readCount ) {
        ArtificialResourcePool gen = new ArtificialResourcePool(createArtificialSamHeader(STARTING_CHROMO, ENDING_CHROMO, READ_COUNT, UNMAPPED_READ_COUNT),
                ArtificialSAMUtils.mappedAndUnmappedReadIterator(STARTING_CHROMO, ENDING_CHROMO, READ_COUNT, UNMAPPED_READ_COUNT));

        GenomeLocParser.setupRefContigOrdering(gen.getHeader().getSequenceDictionary());
        int unmappedReadsSeen = 0;
        int iterations = 0;

        IndexDrivenSAMDataSource data = new IndexDrivenSAMDataSource(reads);
        data.setResourcePool(gen);
        GenomeLocSortedSet set = new GenomeLocSortedSet();
        set.add(GenomeLocParser.createGenomeLoc(0, start, stop));
        ShardStrategy strat = ShardStrategyFactory.shatter(data,null,ShardStrategyFactory.SHATTER_STRATEGY.INTERVAL, gen.getHeader().getSequenceDictionary(), UNMAPPED_READ_COUNT, set);

        StingSAMIterator iter = data.seek(strat.next());
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord r = iter.next();
            // uncomment for debugging - System.err.println(r.getAlignmentStart() + " " + r.getAlignmentEnd());
            count++;
        }
        assertEquals(readCount, count);

    }

    /**
     * test out that we get a single read, given the specific size
     */
    @Test
    public void testSingleRead() {
        testRead(1,ArtificialSAMUtils.DEFAULT_READ_LENGTH,50);
    }

    /**
     * test out that we get the expected amount for a whole chromosome
     */
    @Test
    public void testChromosome() {
        testRead(1, READ_COUNT, READ_COUNT);
    }

    /**
     * test out that we get the expected amount for a whole chromosome
     */
    @Test
    public void testMiddle() {
        testRead(20, READ_COUNT-20,61);
    }


    private SAMFileHeader createArtificialSamHeader( int startingChr, int endingChr, int readCount, int readSize ) {
        return ArtificialSAMUtils.createArtificialSamHeader(( endingChr - startingChr ) + 1,
                startingChr,
                readCount + readSize);
    }
}
