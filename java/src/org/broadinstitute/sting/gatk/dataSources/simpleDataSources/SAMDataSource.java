package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.dataSources.shards.ReadShard;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 2:36:16 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class SAMDataSource implements SimpleDataSource {


    /** Backing support for reads. */
    private Reads reads = null;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMDataSource.class);

    // used for the reads case, the last count of reads retrieved
    long readsTaken = 0;

    // our last genome loc position
    GenomeLoc lastReadPos = null;

    // do we take unmapped reads
    private boolean includeUnmappedReads = true;

    // reads based traversal variables
    private boolean intoUnmappedReads = false;
    private int readsSeenAtLastPos = 0;

    // A pool of SAM iterators.
    private SAMIteratorPool iteratorPool = null;

    /**
     * constructor, given sam files
     *
     * @param reads the list of sam files
     * @param byReads are we a by reads traversal, or a loci traversal.  We could delete this field
     *                if we passed in iterGen, which would be a better (although more complicated for the
     *                consumers of SAMDataSources).
     */
    public SAMDataSource( Reads reads, boolean byReads ) throws SimpleDataSourceLoadException {
        this.reads = reads;

        // check the length
        if (reads.getReadsFiles().size() < 1) {
            throw new SimpleDataSourceLoadException("SAMDataSource: you must provide a list of length greater then 0");
        }
        for (File smFile : reads.getReadsFiles()) {
            if (!smFile.canRead()) {
                throw new SimpleDataSourceLoadException("SAMDataSource: Unable to load file: " + smFile.getName());
            }
        }
        iteratorPool = new SAMIteratorPool(reads,byReads);
    }

    /**
     * For unit testing, add a custom iterator pool.
     * @param iteratorPool Custom mock iterator pool.
     */
    void setResourcePool( SAMIteratorPool iteratorPool ) {
        this.iteratorPool = iteratorPool;
    }


    /**
     * <p>
     * seekLocus
     * </p>
     *
     * @param location the genome location to extract data for
     *
     * @return an iterator for that region
     */
    public StingSAMIterator seekLocus(GenomeLoc location) throws SimpleDataSourceLoadException {
        return iteratorPool.iterator(location);
    }

    /**
     * <p>
     * seek
     * </p>
     *
     * @param shard the shard to get data for
     *
     * @return an iterator for that region
     */
    public StingSAMIterator seek( Shard shard ) throws SimpleDataSourceLoadException {
        StingSAMIterator iterator = null;
        if (shard.getShardType() == Shard.ShardType.READ) {
            iterator = seekRead((ReadShard) shard);
            iterator = TraversalEngine.applyDecoratingIterators(true,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getMaxOnTheFlySorts(),
                    reads.getFilterZeroMappingQualityReads(),
                    reads.getSafetyChecking());
        } else if (shard.getShardType() == Shard.ShardType.LOCUS || shard.getShardType() == Shard.ShardType.INTERVAL) {
            iterator = seekLocus(shard.getGenomeLoc());
            iterator = TraversalEngine.applyDecoratingIterators(false,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getMaxOnTheFlySorts(),
                    reads.getFilterZeroMappingQualityReads(),
                    reads.getSafetyChecking());
        } else {
            throw new StingException("seek: Unknown shard type");
        }

        return iterator;
    }


    /**
     * Gets the (potentially merged) SAM file header.
     *
     * @return SAM file header.
     */
    public SAMFileHeader getHeader() {
        return iteratorPool.getHeader();
    }


    /**
     * <p>
     * seek
     * </p>
     *
     * @param shard the read shard to extract from
     *
     * @return an iterator for that region
     */
    private BoundedReadIterator seekRead( ReadShard shard ) throws SimpleDataSourceLoadException {

        BoundedReadIterator bound = null;
        StingSAMIterator iter = null;

        if (!intoUnmappedReads) {
            if (lastReadPos == null) {
                lastReadPos = new GenomeLoc(getHeader().getSequenceDictionary().getSequence(0).getSequenceIndex(), 0, Integer.MAX_VALUE);
                iter = iteratorPool.iterator(lastReadPos);
                return InitialReadIterator(shard.getSize(), iter);
            } else {
                lastReadPos.setStop(-1);
                iter = iteratorPool.iterator(lastReadPos);
                bound = fastMappedReadSeek(shard.getSize(), StingSAMIteratorAdapter.adapt(reads, iter));
            }
        }

        if (( bound == null || intoUnmappedReads ) && includeUnmappedReads) {
            if (iter != null) {
                iter.close();
            }
            iter = iteratorPool.iterator(null);
            bound = toUnmappedReads(shard.getSize(), iter);
        }
        if (bound == null) {
            shard.signalDone();
            bound = new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), 0);
        }
        return bound;
    }

    /**
     * If we're in by-read mode, this indicates if we want
     * to see unmapped reads too.  Only seeing mapped reads
     * is much faster, but most BAM files have significant
     * unmapped read counts.
     *
     * @param seeUnMappedReads true to see unmapped reads, false otherwise
     */
    public void viewUnmappedReads( boolean seeUnMappedReads ) {
        includeUnmappedReads = seeUnMappedReads;
    }

    /**
     * Seek, if we want unmapped reads.  This method will be faster then the unmapped read method, but you cannot extract the
     * unmapped reads.
     *
     * @param readCount how many reads to retrieve
     * @param iter      the iterator to use
     *
     * @return the bounded iterator that you can use to get the intervaled reads from
     * @throws SimpleDataSourceLoadException
     */
    BoundedReadIterator toUnmappedReads( long readCount, StingSAMIterator iter ) throws SimpleDataSourceLoadException {
        PeekableIterator<SAMRecord> peekable = new PeekableIterator<SAMRecord>(iter);

        int count = 0;
        int cnt = 0;
        SAMRecord d = null;
        while (peekable.hasNext()) {
            d = peekable.peek();
            int x = d.getReferenceIndex();
            if (x < 0)
                // we have the magic read that starts the unmapped read segment!
                break;
            cnt++;
            peekable.next();
        }

        // check to see what happened, did we run out of reads?
        if (!peekable.hasNext()) {
            return null;
        }

        // now walk until we've taken the unmapped read count
        while (peekable.hasNext() && count < this.readsTaken) {
            peekable.next();
            count++;
        }

        // check to see what happened, did we run out of reads?
        if (!peekable.hasNext()) {
            return null;
        }

        // we're not out of unmapped reads, so increment our read cout
        this.readsTaken += readCount;
        return new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, peekable), readCount);

    }


    /**
     * A seek function for unmapped reads.
     *
     * @param readCount how many reads to retrieve
     * @param iter      the iterator to use, seeked to the correct start location
     *
     * @return the bounded iterator that you can use to get the intervaled reads from
     * @throws SimpleDataSourceLoadException
     */
    BoundedReadIterator fastMappedReadSeek( long readCount, StingSAMIterator iter ) throws SimpleDataSourceLoadException {
        BoundedReadIterator bound;
        correctForReadPileupSeek(iter);
        if (readsTaken == 0) {
            return InitialReadIterator(readCount, iter);
        }
        int x = 0;
        SAMRecord rec = null;
        int lastPos = 0;

        while (x < readsTaken) {
            if (iter.hasNext()) {
                rec = iter.next();
                if (lastPos == rec.getAlignmentStart()) ++this.readsSeenAtLastPos;
                else this.readsSeenAtLastPos = 1;
                lastPos = rec.getAlignmentStart();
                ++x;
            } else {
                // jump contigs
                if (lastReadPos.toNextContig() == false) {
                    // check to see if we're using unmapped reads, if not return, we're done
                    readsTaken = 0;
                    intoUnmappedReads = true;
                    return null;
                } else {
                    readsTaken = readCount;
                    readsSeenAtLastPos = 0;
                    lastReadPos.setStop(-1);
                    CloseableIterator<SAMRecord> ret = iteratorPool.iterator(lastReadPos);
                    return new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, ret), readCount);
                }
            }
        }

        // if we're off the end of the last contig (into unmapped territory)
        if (rec != null && rec.getAlignmentStart() == 0) {
            readsTaken += readCount;
            intoUnmappedReads = true;
        }
        // else we're not off the end, store our location
        else if (rec != null) {
            int stopPos = rec.getAlignmentStart();
            if (stopPos < lastReadPos.getStart()) {
                lastReadPos = new GenomeLoc(lastReadPos.getContigIndex() + 1, stopPos, stopPos);
            } else {
                lastReadPos.setStart(rec.getAlignmentStart());
            }
        }
        // in case we're run out of reads, get out
        else {
            throw new StingException("Danger: weve run out reads in fastMappedReadSeek");
            //return null;
        }
        bound = new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), readCount);


        // return the iterator
        return bound;
    }

    /**
     * Even though the iterator has seeked to the correct location, there may be multiple reads at that location,
     * and we may have given some of them out already.  Move the iterator to the correct location using the readsAtLastPos variable
     *
     * @param iter the iterator
     */
    private void correctForReadPileupSeek( StingSAMIterator iter ) {
        // move the number of reads we read from the last pos
        boolean atLeastOneReadSeen = false; // we have a problem where some chomesomes don't have a single read (i.e. the chrN_random chrom.)
        while (iter.hasNext() && this.readsSeenAtLastPos > 0) {
            iter.next();
            --readsSeenAtLastPos;
            atLeastOneReadSeen = true;
        }
        if (readsSeenAtLastPos != 0 && atLeastOneReadSeen) {
            throw new SimpleDataSourceLoadException("Seek problem: reads at last position count != 0");
        }
    }


    /**
     * set the initial iterator
     *
     * @param readCount the number of reads
     * @param iter      the merging iterator
     *
     * @return a bounded read iterator at the first read position in the file.
     */
    private BoundedReadIterator InitialReadIterator( long readCount, CloseableIterator<SAMRecord> iter ) {
        BoundedReadIterator bound;
        bound = new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), readCount);
        this.readsTaken = readCount;
        return bound;
    }

}

class SAMIteratorPool extends ResourcePool<SamFileHeaderMerger,StingSAMIterator> {
    /**
     * Source information about the reads.
     */
    protected Reads reads;

    /**
     * Is this a by-reads traversal or a by-locus?
     */
    protected boolean byReads;

    /**
     * File header for the combined file.
     */
    protected SAMFileHeader header;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMIteratorPool.class);    

    public SAMIteratorPool( Reads reads, boolean byReads ) {
        this.reads = reads;
        this.byReads = byReads;

        SamFileHeaderMerger merger = createNewResource( null );
        this.header = merger.getMergedHeader();
        // Add this resource to the pool.
        this.addNewResource( merger );
    }

    /**
     * Get the combined header for all files in the iterator pool.
     */
    public SAMFileHeader getHeader() {
        return header;
    }

    protected SamFileHeaderMerger selectBestExistingResource( GenomeLoc position, List<SamFileHeaderMerger> mergers) {
        if( mergers.size() == 0 )
            return null;
        return mergers.get(0);
    }

    protected SamFileHeaderMerger createNewResource( GenomeLoc position ) {
        return createHeaderMerger( reads, SAMFileHeader.SortOrder.coordinate );
    }

    protected StingSAMIterator createIteratorFromResource( GenomeLoc loc, SamFileHeaderMerger headerMerger ) {
        final MergingSamRecordIterator2 iterator = new MergingSamRecordIterator2(headerMerger, reads);

        if( loc != null ) {
            if (byReads)
                iterator.queryContained(loc.getContig(), (int) loc.getStart(), (int) loc.getStop());
            else
                iterator.queryOverlapping(loc.getContig(), (int) loc.getStart(), (int) loc.getStop());
        }

        return new StingSAMIterator() {
            public Reads getSourceInfo() { return reads; }
            public void close() {
                iterator.close();
                release(this);
            }
            public Iterator<SAMRecord> iterator() { return this; }
            public boolean hasNext() { return iterator.hasNext(); }
            public SAMRecord next() { return iterator.next(); }
            public void remove() { throw new UnsupportedOperationException("Can't remove from a StingSAMIterator"); }
        };
    }

    protected void closeResource( SamFileHeaderMerger resource ) {
        for( SAMFileReader reader: resource.getReaders() )
            reader.close();
    }

    /**
     * Load a SAM/BAM, given an input file.
     *
     * @param samFile the file name
     *
     * @return a SAMFileReader for the file, null if we're attempting to read a list
     */
    protected SAMFileReader initializeSAMFile( final File samFile, SAMFileReader.ValidationStringency strictness ) {
        if (samFile.toString().endsWith(".list")) {
            return null;
        } else {
            SAMFileReader samReader = new SAMFileReader(samFile, true);
            samReader.setValidationStringency(strictness);

            final SAMFileHeader header = samReader.getFileHeader();
            logger.debug(String.format("Sort order is: " + header.getSortOrder()));

            return samReader;
        }
    }

    /**
     * A private function that, given the internal file list, generates a SamFileReader
     * list of validated files.
     *
     * @return a list of SAMFileReaders that represent the stored file names
     * @throws SimpleDataSourceLoadException if there's a problem loading the files
     */
    protected List<SAMFileReader> GetReaderList( Reads reads ) throws SimpleDataSourceLoadException {
        // right now this is pretty damn heavy, it copies the file list into a reader list every time
        List<SAMFileReader> lst = new ArrayList<SAMFileReader>();
        for (File f : reads.getReadsFiles()) {
            SAMFileReader reader = initializeSAMFile(f, reads.getValidationStringency());

            if (reader.getFileHeader().getReadGroups().size() < 1) {
                //logger.warn("Setting header in reader " + f.getName());
                SAMReadGroupRecord rec = new SAMReadGroupRecord(f.getName());
                rec.setLibrary(f.getName());
                rec.setSample(f.getName());

                reader.getFileHeader().addReadGroup(rec);
            }

            if (reader == null) {
                throw new SimpleDataSourceLoadException("SAMDataSource: Unable to load file: " + f);
            }
            lst.add(reader);
        }
        return lst;
    }

    /**
     * create the merging header.
     *
     * @return a SamFileHeaderMerger that includes the set of SAM files we were created with
     */
    protected SamFileHeaderMerger createHeaderMerger( Reads reads, SAMFileHeader.SortOrder SORT_ORDER ) {
        List<SAMFileReader> lst = GetReaderList(reads);
        return new SamFileHeaderMerger(lst, SORT_ORDER, true);
    }
}