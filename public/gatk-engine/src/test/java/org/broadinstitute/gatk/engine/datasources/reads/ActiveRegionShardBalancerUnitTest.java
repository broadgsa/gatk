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

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileSpan;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;
import java.util.*;

public class ActiveRegionShardBalancerUnitTest extends BaseTest {
    // example genome loc parser for this test, can be deleted if you don't use the reference
    private GenomeLocParser genomeLocParser;
    protected SAMDataSource readsDataSource;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        // sequence
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(10, 0, 10000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        readsDataSource = null;
    }

    @Test
    public void testMergingManyContigs() {
        executeTest(genomeLocParser.getContigs().getSequences());
    }

    @Test
    public void testMergingAllPointersOnSingleContig() {
        executeTest(Arrays.asList(genomeLocParser.getContigs().getSequences().get(1)));
    }

    @Test
    public void testMergingMultipleDiscontinuousContigs() {
        final List<SAMSequenceRecord> all = genomeLocParser.getContigs().getSequences();
        executeTest(Arrays.asList(all.get(1), all.get(3)));
    }

    private void executeTest(final Collection<SAMSequenceRecord> records) {
        final ActiveRegionShardBalancer balancer = new ActiveRegionShardBalancer();

        final List<Set<GenomeLoc>> expectedLocs = new LinkedList<>();
        final List<FilePointer> pointers = new LinkedList<>();

        for ( final SAMSequenceRecord record : records ) {
            final int size = 10;
            int end = 0;
            for ( int i = 0; i < record.getSequenceLength(); i += size) {
                final int myEnd = i + size - 1;
                end = myEnd;
                final GenomeLoc loc = genomeLocParser.createGenomeLoc(record.getSequenceName(), i, myEnd);
                final Map<SAMReaderID, SAMFileSpan> fileSpans = Collections.emptyMap();
                final FilePointer fp = new FilePointer(fileSpans, IntervalMergingRule.ALL, Collections.singletonList(loc));
                pointers.add(fp);
            }
            expectedLocs.add(Collections.singleton(genomeLocParser.createGenomeLoc(record.getSequenceName(), 0, end)));
        }

        balancer.initialize(readsDataSource, pointers.iterator(), genomeLocParser);

        int i = 0;
        int nShardsFound = 0;
        for ( final Shard shard : balancer ) {
            nShardsFound++;
            Assert.assertEquals(new HashSet<>(shard.getGenomeLocs()), expectedLocs.get(i++));
        }
        Assert.assertEquals(nShardsFound, records.size(), "Didn't find exactly one shard for each contig in the sequence dictionary");
    }
}
