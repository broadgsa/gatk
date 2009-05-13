package org.broadinstitute.sting.gatk.dataSources.providers;

import org.junit.Test;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.LocusShard;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.BaseTest;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileHeader;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Collections;
import java.io.FileNotFoundException;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
/**
 * User: hanna
 * Date: May 12, 2009
 * Time: 2:34:46 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Test the locus context queue.
 */
public class LocusContextQueueTest extends BaseTest {
    private static ReferenceSequenceFile sequenceSourceFile = null;

    @BeforeClass
    public static void setupGenomeLoc() throws FileNotFoundException {
        sequenceSourceFile = fakeReferenceSequenceFile();
        GenomeLoc.setupRefContigOrdering(sequenceSourceFile);
    }

    private static ReferenceSequenceFile fakeReferenceSequenceFile() {
        return new ReferenceSequenceFile() {
            public SAMSequenceDictionary getSequenceDictionary() {
                SAMSequenceRecord sequenceRecord = new SAMSequenceRecord("chr1");
                SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Collections.singletonList(sequenceRecord));
                return dictionary;
            }
            public ReferenceSequence nextSequence() { throw new UnsupportedOperationException("Fake implementation doesn't support a getter"); }
        };
    }

    @Test
    public void emptyLocusContextTest() {
        SAMRecordIterator iterator = new SAMRecordIterator();

        GenomeLoc shardBounds = new GenomeLoc("chr1",1,5);
        Shard shard = new LocusShard(shardBounds);
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );

        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        testReadsInContext( queue, shard.getGenomeLoc(), Collections.<SAMRecord>emptyList() );
    }

    @Test
    public void singleReadTest() {
        SAMRecord read = buildSAMRecord("chr1",1,5);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        GenomeLoc shardBounds = new GenomeLoc("chr1",1,5);
        Shard shard = new LocusShard(shardBounds);
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );

        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        testReadsInContext( queue, shard.getGenomeLoc(), Collections.singletonList(read) );        
    }

    @Test
    public void readCoveringFirstPartTest() {
        SAMRecord read = buildSAMRecord("chr1",1,5);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        testReadsInContext( queue, shard.getGenomeLoc(), Collections.singletonList(read) );
    }

    @Test
    public void readCoveringLastPartTest() {
        SAMRecord read = buildSAMRecord("chr1",6,10);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        testReadsInContext( queue, shard.getGenomeLoc(), Collections.singletonList(read) );
    }

    @Test
    public void readCoveringMiddleTest() {
        SAMRecord read = buildSAMRecord("chr1",3,7);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        testReadsInContext( queue, shard.getGenomeLoc(), Collections.singletonList(read) );
    }

    @Test
    public void readOverlappingStartTest() {
        SAMRecord read = buildSAMRecord("chr1",1,10);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new LocusShard(new GenomeLoc("chr1",6,15));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        testReadsInContext( queue, shard.getGenomeLoc(), Collections.singletonList(read) );
    }

    @Test
    public void readOverlappingEndTest() {
        SAMRecord read = buildSAMRecord("chr1",6,15);
        SAMRecordIterator iterator = new SAMRecordIterator(read);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        testReadsInContext( queue, shard.getGenomeLoc(), Collections.singletonList(read) );
    }

    @Test
    public void readsSpanningTest() {
        SAMRecord read1 = buildSAMRecord("chr1",1,5);
        SAMRecord read2 = buildSAMRecord("chr1",6,10);
        SAMRecordIterator iterator = new SAMRecordIterator(read1,read2);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        List<SAMRecord> expectedReads = new ArrayList<SAMRecord>();
        Collections.addAll(expectedReads,read1,read2);
        testReadsInContext( queue, shard.getGenomeLoc(), expectedReads );
    }

    @Test
    public void duplicateReadsTest() {
        SAMRecord read1 = buildSAMRecord("chr1",1,5);
        SAMRecord read2 = buildSAMRecord("chr1",1,5);
        SAMRecord read3 = buildSAMRecord("chr1",6,10);
        SAMRecord read4 = buildSAMRecord("chr1",6,10);
        SAMRecordIterator iterator = new SAMRecordIterator(read1,read2,read3,read4);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        List<SAMRecord> expectedReads = new ArrayList<SAMRecord>();
        Collections.addAll(expectedReads,read1,read2,read3,read4);
        testReadsInContext( queue, shard.getGenomeLoc(), expectedReads );
    }

    @Test
    public void cascadingReadsWithinBoundsTest() {
        SAMRecord read1 = buildSAMRecord("chr1",2,6);
        SAMRecord read2 = buildSAMRecord("chr1",3,7);
        SAMRecord read3 = buildSAMRecord("chr1",4,8);
        SAMRecord read4 = buildSAMRecord("chr1",5,9);
        SAMRecordIterator iterator = new SAMRecordIterator(read1,read2,read3,read4);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        List<SAMRecord> expectedReads = new ArrayList<SAMRecord>();
        Collections.addAll(expectedReads,read1,read2,read3,read4);
        testReadsInContext( queue, shard.getGenomeLoc(), expectedReads );
    }

    @Test
    public void cascadingReadsAtBoundsTest() {
        SAMRecord read1 = buildSAMRecord("chr1",1,5);
        SAMRecord read2 = buildSAMRecord("chr1",2,6);
        SAMRecord read3 = buildSAMRecord("chr1",3,7);
        SAMRecord read4 = buildSAMRecord("chr1",4,8);
        SAMRecord read5 = buildSAMRecord("chr1",5,9);
        SAMRecord read6 = buildSAMRecord("chr1",6,10);
        SAMRecordIterator iterator = new SAMRecordIterator(read1,read2,read3,read4,read5,read6);

        Shard shard = new LocusShard(new GenomeLoc("chr1",1,10));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        List<SAMRecord> expectedReads = new ArrayList<SAMRecord>();
        Collections.addAll(expectedReads,read1,read2,read3,read4,read5,read6);
        testReadsInContext( queue, shard.getGenomeLoc(), expectedReads );
    }

    @Test
    public void cascadingReadsOverlappingBoundsTest() {
        SAMRecord read01 = buildSAMRecord("chr1",1,5);
        SAMRecord read02 = buildSAMRecord("chr1",2,6);
        SAMRecord read03 = buildSAMRecord("chr1",3,7);
        SAMRecord read04 = buildSAMRecord("chr1",4,8);
        SAMRecord read05 = buildSAMRecord("chr1",5,9);
        SAMRecord read06 = buildSAMRecord("chr1",6,10);
        SAMRecord read07 = buildSAMRecord("chr1",7,11);
        SAMRecord read08 = buildSAMRecord("chr1",8,12);
        SAMRecord read09 = buildSAMRecord("chr1",9,13);
        SAMRecord read10 = buildSAMRecord("chr1",10,14);
        SAMRecord read11 = buildSAMRecord("chr1",11,15);
        SAMRecord read12 = buildSAMRecord("chr1",12,16);
        SAMRecordIterator iterator = new SAMRecordIterator(read01,read02,read03,read04,read05,read06,
                                                           read07,read08,read09,read10,read11,read12);

        Shard shard = new LocusShard(new GenomeLoc("chr1",6,15));
        ShardDataProvider dataProvider = new ShardDataProvider( shard, iterator );
        LocusContextQueue queue = new LocusContextQueue( dataProvider );

        List<SAMRecord> expectedReads = new ArrayList<SAMRecord>();
        Collections.addAll(expectedReads,read01,read02,read03,read04,read05,read06,
                                         read07,read08,read09,read10,read11,read12);
        testReadsInContext( queue, shard.getGenomeLoc(), expectedReads );
    }

    /**
     * Test the reads according to an independently derived context.
     * @param queue
     * @param bounds
     * @param reads
     */
    private void testReadsInContext( LocusContextQueue queue, GenomeLoc bounds, List<SAMRecord> reads ) {
        Assert.assertEquals("Initial position of queue is incorrect", new GenomeLoc(bounds.getContig(),bounds.getStart()), queue.getSeekPoint() );

        for( long i = bounds.getStart(); i <= bounds.getStop(); i++ ) {
            GenomeLoc site = new GenomeLoc("chr1",i);
            queue.seek(site);
            Assert.assertEquals("Seeked queue is incorrect", site, queue.getSeekPoint() );

            LocusContext locusContext = queue.peek();
            Assert.assertEquals("Target locus context location is incorrect", site, locusContext.getLocation() );
            int expectedReadsAtSite = 0;

            for( SAMRecord read: reads ) {
                if(new GenomeLoc(read).containsP(locusContext.getLocation())) {
                    Assert.assertTrue("Target locus context does not contain reads", locusContext.getReads().contains(read) );
                    expectedReadsAtSite++;
                }
            }

            Assert.assertEquals("Found wrong number of reads at site", expectedReadsAtSite, locusContext.getReads().size());
        }

    }

    private SAMRecord buildSAMRecord( String contig, int alignmentStart, int alignmentEnd ) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(sequenceSourceFile.getSequenceDictionary());

        SAMRecord record = new SAMRecord(header);        

        record.setReferenceIndex(sequenceSourceFile.getSequenceDictionary().getSequenceIndex(contig));
        record.setAlignmentStart(alignmentStart);
        Cigar cigar = new Cigar();
        cigar.add(new CigarElement(alignmentEnd-alignmentStart+1,CigarOperator.M));
        record.setCigar(cigar);
        return record;
    }

    private class SAMRecordIterator implements StingSAMIterator {
        private Iterator<SAMRecord> backingIterator = null;

        public SAMRecordIterator( SAMRecord... reads ) {
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
