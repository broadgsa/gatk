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

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.engine.datasources.reads.MockLocusShard;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.engine.executive.WindowMaker;
import org.broadinstitute.gatk.engine.datasources.reads.LocusShard;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
/**
 * User: hanna
 * Date: May 13, 2009
 * Time: 4:29:08 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/** Base support for testing variants of the LocusView family of classes. */

public abstract class LocusViewTemplate extends BaseTest {
    protected static ReferenceSequenceFile sequenceSourceFile = null;
    protected GenomeLocParser genomeLocParser = null;

    @BeforeClass
    public void setupGenomeLoc() throws FileNotFoundException {
        sequenceSourceFile = fakeReferenceSequenceFile();
        genomeLocParser = new GenomeLocParser(sequenceSourceFile);
    }

    @Test
    public void emptyAlignmentContextTest() {
        SAMRecordIterator iterator = new SAMRecordIterator();

        GenomeLoc shardBounds = genomeLocParser.createGenomeLoc("chr1", 1, 5);
        Shard shard = new LocusShard(genomeLocParser, new SAMDataSource(null,Collections.<SAMReaderID>emptyList(),new ThreadAllocation(),null,genomeLocParser),Collections.singletonList(shardBounds),Collections.<SAMReaderID,SAMFileSpan>emptyMap());
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, null, genomeLocParser, window.getLocus(), window, null, null);

        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.<GATKSAMRecord>emptyList());
    }

    @Test
    public void singleReadTest() {
        GATKSAMRecord read = buildSAMRecord("read1","chr1", 1, 5);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        GenomeLoc shardBounds = genomeLocParser.createGenomeLoc("chr1", 1, 5);
        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(shardBounds));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);

        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.singletonList(read));
    }

    @Test
    public void readCoveringFirstPartTest() {
        GATKSAMRecord read = buildSAMRecord("read1","chr1", 1, 5);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.singletonList(read));
    }

    @Test
    public void readCoveringLastPartTest() {
        GATKSAMRecord read = buildSAMRecord("read1","chr1", 6, 10);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.singletonList(read));
    }

    @Test
    public void readCoveringMiddleTest() {
        GATKSAMRecord read = buildSAMRecord("read1","chr1", 3, 7);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.singletonList(read));
    }

    @Test
    public void readAndLocusOverlapAtLastBase() {
        GATKSAMRecord read = buildSAMRecord("read1","chr1", 1, 5);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 5, 5)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.singletonList(read));
    }

    @Test
    public void readOverlappingStartTest() {
        GATKSAMRecord read = buildSAMRecord("read1","chr1", 1, 10);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 6, 15)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.singletonList(read));
    }

    @Test
    public void readOverlappingEndTest() {
        GATKSAMRecord read = buildSAMRecord("read1","chr1", 6, 15);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        testReadsInContext(view, shard.getGenomeLocs(), Collections.singletonList(read));
    }

    @Test
    public void readsSpanningTest() {
        GATKSAMRecord read1 = buildSAMRecord("read1","chr1", 1, 5);
        GATKSAMRecord read2 = buildSAMRecord("read2","chr1", 6, 10);
        SAMRecordIterator iterator = new SAMRecordIterator(read1, read2);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        List<GATKSAMRecord> expectedReads = new ArrayList<GATKSAMRecord>();
        Collections.addAll(expectedReads, read1, read2);
        testReadsInContext(view, shard.getGenomeLocs(), expectedReads);
    }

    @Test
    public void duplicateReadsTest() {
        GATKSAMRecord read1 = buildSAMRecord("read1","chr1", 1, 5);
        GATKSAMRecord read2 = buildSAMRecord("read2","chr1", 1, 5);
        GATKSAMRecord read3 = buildSAMRecord("read3","chr1", 6, 10);
        GATKSAMRecord read4 = buildSAMRecord("read4","chr1", 6, 10);
        SAMRecordIterator iterator = new SAMRecordIterator(read1, read2, read3, read4);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        List<GATKSAMRecord> expectedReads = new ArrayList<GATKSAMRecord>();
        Collections.addAll(expectedReads, read1, read2, read3, read4);
        testReadsInContext(view, shard.getGenomeLocs(), expectedReads);
    }

    @Test
    public void cascadingReadsWithinBoundsTest() {
        GATKSAMRecord read1 = buildSAMRecord("read1","chr1", 2, 6);
        GATKSAMRecord read2 = buildSAMRecord("read2","chr1", 3, 7);
        GATKSAMRecord read3 = buildSAMRecord("read3","chr1", 4, 8);
        GATKSAMRecord read4 = buildSAMRecord("read4","chr1", 5, 9);
        SAMRecordIterator iterator = new SAMRecordIterator(read1, read2, read3, read4);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        List<GATKSAMRecord> expectedReads = new ArrayList<GATKSAMRecord>();
        Collections.addAll(expectedReads, read1, read2, read3, read4);
        testReadsInContext(view, shard.getGenomeLocs(), expectedReads);
    }

    @Test
    public void cascadingReadsAtBoundsTest() {
        GATKSAMRecord read1 = buildSAMRecord("read1","chr1", 1, 5);
        GATKSAMRecord read2 = buildSAMRecord("read2","chr1", 2, 6);
        GATKSAMRecord read3 = buildSAMRecord("read3","chr1", 3, 7);
        GATKSAMRecord read4 = buildSAMRecord("read4","chr1", 4, 8);
        GATKSAMRecord read5 = buildSAMRecord("read5","chr1", 5, 9);
        GATKSAMRecord read6 = buildSAMRecord("read6","chr1", 6, 10);
        SAMRecordIterator iterator = new SAMRecordIterator(read1, read2, read3, read4, read5, read6);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 1, 10)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        List<GATKSAMRecord> expectedReads = new ArrayList<GATKSAMRecord>();
        Collections.addAll(expectedReads, read1, read2, read3, read4, read5, read6);
        testReadsInContext(view, shard.getGenomeLocs(), expectedReads);
    }

    @Test
    public void cascadingReadsOverlappingBoundsTest() {
        GATKSAMRecord read01 = buildSAMRecord("read1","chr1", 1, 5);
        GATKSAMRecord read02 = buildSAMRecord("read2","chr1", 2, 6);
        GATKSAMRecord read03 = buildSAMRecord("read3","chr1", 3, 7);
        GATKSAMRecord read04 = buildSAMRecord("read4","chr1", 4, 8);
        GATKSAMRecord read05 = buildSAMRecord("read5","chr1", 5, 9);
        GATKSAMRecord read06 = buildSAMRecord("read6","chr1", 6, 10);
        GATKSAMRecord read07 = buildSAMRecord("read7","chr1", 7, 11);
        GATKSAMRecord read08 = buildSAMRecord("read8","chr1", 8, 12);
        GATKSAMRecord read09 = buildSAMRecord("read9","chr1", 9, 13);
        GATKSAMRecord read10 = buildSAMRecord("read10","chr1", 10, 14);
        GATKSAMRecord read11 = buildSAMRecord("read11","chr1", 11, 15);
        GATKSAMRecord read12 = buildSAMRecord("read12","chr1", 12, 16);
        SAMRecordIterator iterator = new SAMRecordIterator(read01, read02, read03, read04, read05, read06,
                                                           read07, read08, read09, read10, read11, read12);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chr1", 6, 15)));
        WindowMaker windowMaker = new WindowMaker(shard,genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();
        LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, window.getSourceInfo(), genomeLocParser, window.getLocus(), window, null, null);
        LocusView view = createView(dataProvider);

        List<GATKSAMRecord> expectedReads = new ArrayList<GATKSAMRecord>();
        Collections.addAll(expectedReads, read01, read02, read03, read04, read05, read06,
                           read07, read08, read09, read10, read11, read12);
        testReadsInContext(view, shard.getGenomeLocs(), expectedReads);
    }

    /**
     * Creates a view of the type required for testing.
     *
     * @return The correct view to test.
     */
    protected abstract LocusView createView(LocusShardDataProvider provider);

    /**
     * Test the reads according to an independently derived context.
     *
     * @param view
     * @param bounds
     * @param reads
     */
    protected abstract void testReadsInContext(LocusView view, List<GenomeLoc> bounds, List<GATKSAMRecord> reads);

    /**
     * Fake a reference sequence file.  Essentially, seek a header with a bunch of dummy data.
     *
     * @return A 'fake' reference sequence file
     */
    private static ReferenceSequenceFile fakeReferenceSequenceFile() {
        return new ReferenceSequenceFile() {
            public SAMSequenceDictionary getSequenceDictionary() {
                SAMSequenceRecord sequenceRecord = new SAMSequenceRecord("chr1", 1000000);
                SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Collections.singletonList(sequenceRecord));
                return dictionary;
            }

            public boolean isIndexed() { return false; }

            public ReferenceSequence nextSequence() {
                throw new UnsupportedOperationException("Fake implementation doesn't support a getter");
            }

            public ReferenceSequence getSequence( String contig ) {
                throw new UnsupportedOperationException("Fake implementation doesn't support a getter");
            }

            public ReferenceSequence getSubsequenceAt( String contig, long start, long stop ) {
                throw new UnsupportedOperationException("Fake implementation doesn't support a getter");
            }

            public void reset() {
                return;
            }

            public void close() throws IOException {
            }
        };
    }

    /**
     * Build a SAM record featuring the absolute minimum required dataset.
     *
     * @param contig         Contig to populate.
     * @param alignmentStart start of alignment
     * @param alignmentEnd   end of alignment
     *
     * @return New SAM Record
     */
    protected GATKSAMRecord buildSAMRecord(String readName, String contig, int alignmentStart, int alignmentEnd) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(sequenceSourceFile.getSequenceDictionary());

        GATKSAMRecord record = new GATKSAMRecord(header);

        record.setReadName(readName);
        record.setReferenceIndex(sequenceSourceFile.getSequenceDictionary().getSequenceIndex(contig));
        record.setAlignmentStart(alignmentStart);
        Cigar cigar = new Cigar();
        int len = alignmentEnd - alignmentStart + 1;
        cigar.add(new CigarElement(len, CigarOperator.M));
        record.setCigar(cigar);
        record.setReadBases(new byte[len]);
        record.setBaseQualities(new byte[len]);
        return record;
    }

    /** A simple iterator which iterates over a list of reads. */
    protected class SAMRecordIterator implements GATKSAMIterator {
        private Iterator<SAMRecord> backingIterator = null;

        public SAMRecordIterator(SAMRecord... reads) {
            List<SAMRecord> backingList = new ArrayList<SAMRecord>();
            backingList.addAll(Arrays.asList(reads));
            backingIterator = backingList.iterator();
        }

        public boolean hasNext() {
            return backingIterator.hasNext();
        }

        public SAMRecord next() {
            return backingIterator.next();
        }

        public Iterator<SAMRecord> iterator() {
            return this;
        }

        public void close() {
            // NO-OP.
        }

        public void remove() {
            throw new UnsupportedOperationException("Can't remove from a read-only iterator");
        }
    }


}
