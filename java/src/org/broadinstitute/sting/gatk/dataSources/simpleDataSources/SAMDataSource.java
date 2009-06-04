package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import net.sf.picard.sam.SamFileHeaderMerger;
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


    /**
     * package protected getter and setter for the iterator generator
     *
     * @return
     */
    IteratorGenerator getIterGen() {
        return iterGen;
    }

    void setIterGen( IteratorGenerator iterGen ) {
        this.iterGen = iterGen;
    }

    // where we get out iterators from
    private IteratorGenerator iterGen;

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
        iterGen = new MSR2IteratorGenerator(reads, byReads);

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
    public StingSAMIterator seekLocus( GenomeLoc location ) throws SimpleDataSourceLoadException {
       // make a merging iterator for this record
        CloseableIterator<SAMRecord> iter = iterGen.seek(location);
        // return the iterator
        return StingSAMIteratorAdapter.adapt(reads, iter);
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
        return this.iterGen.getHeader();
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
        CloseableIterator<SAMRecord> iter = null;

        if (!intoUnmappedReads) {
            if (lastReadPos == null) {
                lastReadPos = new GenomeLoc(iterGen.getHeader().getSequenceDictionary().getSequence(0).getSequenceIndex(), 0, Integer.MAX_VALUE);
                iter = iterGen.seek(lastReadPos);
                return InitialReadIterator(shard.getSize(), iter);
            } else {
                lastReadPos.setStop(-1);
                iter = iterGen.seek(lastReadPos);
                bound = fastMappedReadSeek(shard.getSize(), StingSAMIteratorAdapter.adapt(reads, iter));
            }
        }

        if (( bound == null || intoUnmappedReads ) && includeUnmappedReads) {
            if (iter != null) {
                iter.close();
            }
            iter = iterGen.seek(null);
            bound = toUnmappedReads(shard.getSize(), (PeekingStingIterator) iter);
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
    <U extends PeekingStingIterator> BoundedReadIterator toUnmappedReads( long readCount, U iter ) throws SimpleDataSourceLoadException {
        int count = 0;
        int cnt = 0;
        SAMRecord d = null;
        while (iter.hasNext()) {
            d = iter.peek();
            int x = d.getReferenceIndex();
            if (x < 0)
                // we have the magic read that starts the unmapped read segment!
                break;
            cnt++;
            iter.next();
        }

        // check to see what happened, did we run out of reads?
        if (!iter.hasNext()) {
            return null;
        }

        // now walk until we've taken the unmapped read count
        while (iter.hasNext() && count < this.readsTaken) {
            iter.next();
            count++;
        }

        // check to see what happened, did we run out of reads?
        if (!iter.hasNext()) {
            return null;
        }

        // we're not out of unmapped reads, so increment our read cout
        this.readsTaken += readCount;
        return new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), readCount);

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
                    CloseableIterator<SAMRecord> ret = iterGen.seek(lastReadPos);
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


/**
 * iterator generator
 * <p/>
 * This class generates iterators for the SAMDataSource.  This class was introduced for testing purposes,
 * since it became increasingly hard to test the SAM data source code.  The class defines two abstraact
 * methods:
 * <p/>
 * -seek( GenomeLoc ) which returns an iterator seeked to the genome loc, and if null is passed to the default
 * location (which is implementation specific).
 * <p/>
 * -getHeader(), which returns a SAMFileHeader for the specified IteratorGenerator.  I hope we can phase this
 * method out, since it doesn't seem necessary, and it would be much cleaner with out it.
 */
abstract class IteratorGenerator {
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMDataSource.class);

    /**
     * seek to a location
     *
     * @param seekTo the genome loc to seek to
     *
     * @return StingSAMIterator
     */
    public abstract CloseableIterator<SAMRecord> seek( GenomeLoc seekTo );

    /**
     * get the merged header
     *
     * @return the merged header
     */
    public abstract SAMFileHeader getHeader();

    /**
     * Load a SAM/BAM, given an input file.
     *
     * @param samFile the file name
     *
     * @return a SAMFileReader for the file, null if we're attempting to read a list
     */
    protected static SAMFileReader initializeSAMFile( final File samFile, SAMFileReader.ValidationStringency strictness ) {
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
    protected static List<SAMFileReader> GetReaderList( Reads reads, SAMFileReader.ValidationStringency strictness ) throws SimpleDataSourceLoadException {
        // right now this is pretty damn heavy, it copies the file list into a reader list every time
        List<SAMFileReader> lst = new ArrayList<SAMFileReader>();
        for (File f : reads.getReadsFiles()) {
            SAMFileReader reader = initializeSAMFile(f, strictness);

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
    protected static SamFileHeaderMerger createHeaderMerger( Reads reads, SAMFileReader.ValidationStringency strictness, SAMFileHeader.SortOrder SORT_ORDER ) {
        List<SAMFileReader> lst = GetReaderList(reads, strictness);
        return new SamFileHeaderMerger(lst, SORT_ORDER, true);
    }
}


/**
 * MSR2IteratorGenerator
 * <p/>
 * generates a MerginsSAMIterator2, given a genomic location.  The constructor takes the reads structure,
 * and a flag indicating if we're dealing with reads or loci (to determine the correct query function).
 */
class MSR2IteratorGenerator extends IteratorGenerator {
    /** our read pile */
    private Reads reads;

    private SamFileHeaderMerger header;

    // How strict should we be with SAM/BAM parsing?
    protected SAMFileReader.ValidationStringency strictness = SAMFileReader.ValidationStringency.SILENT;

    // are we by reads or by loci
    protected boolean byReads = true;

    /** our SAM data files */
    private final SAMFileHeader.SortOrder sortOrder = SAMFileHeader.SortOrder.coordinate;

    public MSR2IteratorGenerator( Reads reads, boolean byReads ) {
        this.reads = reads;
        this.header = IteratorGenerator.createHeaderMerger(reads, strictness, sortOrder);
        this.byReads = byReads;
    }

    public CloseableIterator<SAMRecord> seek( GenomeLoc seekTo ) {
        SamFileHeaderMerger mg = createHeaderMerger(reads, strictness, sortOrder);
        MergingSamRecordIterator2 iter = new MergingSamRecordIterator2(mg, reads);
        if (seekTo != null) {
            if (byReads)
                iter.queryContained(seekTo.getContig(), (int) seekTo.getStart(), (int) seekTo.getStop());
            else
                iter.queryOverlapping(seekTo.getContig(), (int) seekTo.getStart(), (int) seekTo.getStop());
        }
        return iter;
    }

    /**
     * get the merged header
     *
     * @return the merged header
     */
    public SAMFileHeader getHeader() {
        return header.getMergedHeader();
    }
}


