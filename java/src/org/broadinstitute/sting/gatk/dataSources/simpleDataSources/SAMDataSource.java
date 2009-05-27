package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import edu.mit.broad.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.dataSources.shards.ReadShard;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
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

    /** our SAM data files */
    private final SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMDataSource.class);

    // How strict should we be with SAM/BAM parsing?
    protected SAMFileReader.ValidationStringency strictness = SAMFileReader.ValidationStringency.SILENT;

    // our list of readers
    private final List<File> samFileList = new ArrayList<File>();

    /** SAM header file. */
    private final SAMFileHeader header;

    // used for the reads case, the last count of reads retrieved
    private long readsTaken = 0;

    // our last genome loc position
    private GenomeLoc lastReadPos = null;

    // do we take unmapped reads
    private boolean includeUnmappedReads = true;

    // reads based traversal variables
    private boolean intoUnmappedReads = false;
    private int readsAtLastPos = 0;

    /**
     * constructor, given sam files
     *
     * @param reads the list of sam files
     */
    public SAMDataSource(Reads reads) throws SimpleDataSourceLoadException {
        this.reads = reads;

        // check the length
        if (reads.getReadsFiles().size() < 1) {
            throw new SimpleDataSourceLoadException("SAMDataSource: you must provide a list of length greater then 0");
        }
        for (File smFile : reads.getReadsFiles()) {
            if (!smFile.canRead()) {
                throw new SimpleDataSourceLoadException("SAMDataSource: Unable to load file: " + smFile.getName());
            }
            samFileList.add(smFile);

        }

        header = createHeaderMerger().getMergedHeader();
    }

    /**
     * Load a SAM/BAM, given an input file.
     *
     * @param samFile the file name
     * @return a SAMFileReader for the file, null if we're attempting to read a list
     */
    private SAMFileReader initializeSAMFile(final File samFile) {
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
     * <p>
     * seekLocus
     * </p>
     *
     * @param location the genome location to extract data for
     * @return an iterator for that region
     */
    public StingSAMIterator seekLocus(GenomeLoc location) throws SimpleDataSourceLoadException {

        // right now this is very heavy, it copies the file list into a reader list every time
        SamFileHeaderMerger headerMerger = createHeaderMerger();

        // make a merging iterator for this record
        MergingSamRecordIterator2 iter = new MergingSamRecordIterator2(headerMerger);

        iter.queryOverlapping(location.getContig(), (int) location.getStart(), (int) location.getStop() + 1);

        // return the iterator
        return StingSAMIteratorAdapter.adapt(reads, iter);
    }

    /**
     * <p>
     * seek
     * </p>
     *
     * @param shard the shard to get data for
     * @return an iterator for that region
     */
    public StingSAMIterator seek(Shard shard) throws SimpleDataSourceLoadException {
        StingSAMIterator iterator = null;
        if (shard.getShardType() == Shard.ShardType.READ) {
            iterator = seekRead((ReadShard) shard);
            iterator = TraversalEngine.applyDecoratingIterators(true,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getMaxOnTheFlySorts(),
                    reads.getSafetyChecking());
        } else if (shard.getShardType() == Shard.ShardType.LOCUS ||
                   shard.getShardType() == Shard.ShardType.INTERVAL) {
            iterator = seekLocus(shard.getGenomeLoc());
            iterator = TraversalEngine.applyDecoratingIterators(false,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getMaxOnTheFlySorts(),
                    reads.getSafetyChecking());
        }
        else {
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
        return header;
    }

    /**
     * create the merging header.
     *
     * @return a SamFileHeaderMerger that includes the set of SAM files we were created with
     */
    private SamFileHeaderMerger createHeaderMerger() {
        List<SAMFileReader> lst = GetReaderList();
        return new SamFileHeaderMerger(lst, SORT_ORDER);
    }


    /**
     * <p>
     * seek
     * </p>
     *
     * @param shard the read shard to extract from
     * @return an iterator for that region
     */
    private BoundedReadIterator seekRead(ReadShard shard) throws SimpleDataSourceLoadException {

        BoundedReadIterator bound = null;
        SamFileHeaderMerger headerMerger = createHeaderMerger();
        MergingSamRecordIterator2 iter = null;

        if (!intoUnmappedReads) {
            iter = new MergingSamRecordIterator2(headerMerger);
            bound = fastMappedReadSeek(shard.getSize(), iter);
        }
        if ((bound == null || intoUnmappedReads) && includeUnmappedReads) {
            if (iter != null) {
                iter.close();
            }
            iter = new MergingSamRecordIterator2(createHeaderMerger());
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
    public void viewUnmappedReads(boolean seeUnMappedReads) {
        includeUnmappedReads = seeUnMappedReads;
    }

    /**
     * Seek, if we want unmapped reads.  This method will be faster then the unmapped read method, but you cannot extract the
     * unmapped reads.
     *
     * @param readCount how many reads to retrieve
     * @param iter      the iterator to use
     * @return the bounded iterator that you can use to get the intervaled reads from
     * @throws SimpleDataSourceLoadException
     */
    private BoundedReadIterator toUnmappedReads(long readCount, MergingSamRecordIterator2 iter) throws SimpleDataSourceLoadException {
        int count = 0;
        SAMRecord d = null;
        while (iter.hasNext()) {
            d = iter.peek();
            int x = d.getReferenceIndex();
            if (x < 0 || x >= d.getHeader().getSequenceDictionary().getSequences().size()) {
                // we have the magic read that starts the unmapped read segment!
                break;
            }
            iter.next();
        }

        // check to see what happened, did we run out of reads?
        if (!iter.hasNext()) {
            return null;
        }

        // now walk until we've taken the unmapped read count
        while (iter.hasNext() && count < this.readsTaken) {
            iter.next();
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
     * @param iter      the iterator to use
     * @return the bounded iterator that you can use to get the intervaled reads from
     * @throws SimpleDataSourceLoadException
     */
    private BoundedReadIterator fastMappedReadSeek(long readCount, MergingSamRecordIterator2 iter) throws SimpleDataSourceLoadException {
        if (lastReadPos == null) {
            return InitialReadIterator(readCount, iter);
        } else {
            BoundedReadIterator bound;
            iter.queryContained(lastReadPos.getContig(), (int) lastReadPos.getStop(), -1);

            // move the number of reads we read from the last pos
            boolean atLeastOneReadSeen = false; // we have a problem where some chomesomes don't have a single read (i.e. the chrN_random chrom.)
            while (iter.hasNext() && this.readsAtLastPos > 0) {
                iter.next();
                --readsAtLastPos;
                atLeastOneReadSeen = true;
            }
            if (readsAtLastPos != 0 && atLeastOneReadSeen) {
                throw new SimpleDataSourceLoadException("Seek problem: reads at last position count != 0");
            }

            int x = 0;
            SAMRecord rec = null;
            int lastPos = 0;

            while (x < readsTaken) {
                if (iter.hasNext()) {
                    rec = iter.next();
                    if (lastPos == rec.getAlignmentStart()) {
                        ++this.readsAtLastPos;
                    } else {
                        this.readsAtLastPos = 1;
                    }
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
                        iter.close();
                        // now merge the headers
                        // right now this is pretty damn heavy, it copies the file list into a reader list every time
                        SamFileHeaderMerger mg = createHeaderMerger();
                        iter = new MergingSamRecordIterator2(mg);
                        iter.queryContained(lastReadPos.getContig(), 1, Integer.MAX_VALUE);
                        return new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), readCount);
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
                    lastReadPos.setStop(rec.getAlignmentStart());
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


    }

    /**
     * set the initial iterator
     *
     * @param readCount the number of reads
     * @param iter      the merging iterator
     * @return a bounded read iterator at the first read position in the file.
     */
    private BoundedReadIterator InitialReadIterator(long readCount, MergingSamRecordIterator2 iter) {
        BoundedReadIterator bound;
        lastReadPos = new GenomeLoc(iter.getHeader().getSequenceDictionary().getSequence(0).getSequenceIndex(), 0, 0);
        iter.queryContained(lastReadPos.getContig(), 1, -1);
        bound = new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), readCount);
        this.readsTaken = readCount;
        return bound;
    }


    /**
     * A private function that, given the internal file list, generates a SamFileReader
     * list of validated files.
     *
     * @return a list of SAMFileReaders that represent the stored file names
     * @throws SimpleDataSourceLoadException if there's a problem loading the files
     */
    private List<SAMFileReader> GetReaderList() throws SimpleDataSourceLoadException {
        // right now this is pretty damn heavy, it copies the file list into a reader list every time
        List<SAMFileReader> lst = new ArrayList<SAMFileReader>();
        for (File f : this.samFileList) {
            SAMFileReader reader = initializeSAMFile(f);

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

}
