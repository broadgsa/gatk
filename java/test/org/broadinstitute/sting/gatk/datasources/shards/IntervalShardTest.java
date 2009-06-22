package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
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
 *         Class IntervalReadShardTest
 *         <p/>
 *         Tests for the IntervalReadShard class.
 */
public class IntervalShardTest extends BaseTest {

    private IntervalShard intervalShard = null;
    private SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    @Before
    public void setup() {
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }


    @Test
    public void simpleReturn() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, 1, 100);
        intervalShard = new IntervalShard(loc);
        assertTrue(intervalShard.getGenomeLoc().equals(loc));
    }

    @Test
    public void ensureNotReference() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, 1, 100);
        intervalShard = new IntervalShard(loc);
        assertTrue(intervalShard.getGenomeLoc() != loc && intervalShard.getGenomeLoc().equals(loc));
    }

}
