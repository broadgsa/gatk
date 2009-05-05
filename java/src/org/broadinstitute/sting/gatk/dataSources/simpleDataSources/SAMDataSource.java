package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import edu.mit.broad.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.dataSources.shards.ReadShard;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
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
    /** our SAM data files */
    private final SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMDataSource.class);

    // are we set to locus mode or read mode for dividing
    private boolean locusMode = true;

    // How strict should we be with SAM/BAM parsing?
    protected SAMFileReader.ValidationStringency strictness = SAMFileReader.ValidationStringency.SILENT;

    // our list of readers
    private final List<File> samFileList = new ArrayList<File>();

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
     * @param samFiles the list of sam files
     */
    public SAMDataSource(List<?> samFiles) throws SimpleDataSourceLoadException {
        // check the length
        if (samFiles.size() < 1) {
            throw new SimpleDataSourceLoadException("SAMDataSource: you must provide a list of length greater then 0");
        }
        for (Object fileName : samFiles) {
            File smFile;
            if (samFiles.get(0) instanceof String) {
                smFile = new File((String) samFiles.get(0));
            } else if (samFiles.get(0) instanceof File) {
                smFile = (File) fileName;
            } else {
                throw new SimpleDataSourceLoadException("SAMDataSource: unknown samFile list type, must be String or File");
            }

            if (!smFile.canRead()) {
                throw new SimpleDataSourceLoadException("SAMDataSource: Unable to load file: " + fileName);
            }
            samFileList.add(smFile);

        }

    }

    /**
     * Load up a sam file.
     *
     * @param samFile the file name
     * @return a SAMFileReader for the file
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
     * seek
     * </p>
     *
     * @param location the genome location to extract data for
     * @return an iterator for that region
     */
    public MergingSamRecordIterator2 seekLocus(GenomeLoc location) throws SimpleDataSourceLoadException {

        // right now this is pretty damn heavy, it copies the file list into a reader list every time
        SamFileHeaderMerger headerMerger = CreateHeader();

        // make a merging iterator for this record
        MergingSamRecordIterator2 iter = new MergingSamRecordIterator2(headerMerger);


        // we do different things for locus and read modes
        if (locusMode) {
            iter.queryOverlapping(location.getContig(), (int) location.getStart(), (int) location.getStop() + 1);
        } else {
            iter.queryContained(location.getContig(), (int) location.getStart(), (int) location.getStop() + 1);
        }

        // return the iterator
        return iter;
    }

    /**
     * <p>
     * seek
     * </p>
     *
     * @param shard the shard to get data for
     * @return an iterator for that region
     */
    public CloseableIterator<SAMRecord> seek(Shard shard) throws SimpleDataSourceLoadException {
        if (shard.getShardType() == Shard.ShardType.READ) {
            return seekRead((ReadShard) shard);
        } else if (shard.getShardType() == Shard.ShardType.LOCUS) {
            return seekLocus(shard.getGenomeLoc());
        } else {
            throw new StingException("seek: Unknown shard type");
        }
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
     * <p>
     * seek
     * </p>
     *
     * @param shard the read shard to extract from
     * @return an iterator for that region
     */
    private BoundedReadIterator seekRead(ReadShard shard) throws SimpleDataSourceLoadException {

        BoundedReadIterator bound = null;
        SamFileHeaderMerger headerMerger = CreateHeader();
        MergingSamRecordIterator2 iter = null;

        if (!intoUnmappedReads) {
            // make a merging iterator for this record
            iter = new MergingSamRecordIterator2(headerMerger);

           bound = fastMappedReadSeek(shard.getSize(), iter);
        }
        if ((bound == null || intoUnmappedReads) && includeUnmappedReads) {
            if (iter != null) {
                iter.close();
            }
            iter = new MergingSamRecordIterator2(CreateHeader());
            bound = toUnmappedReads(shard.getSize(), iter);
        }

        if (bound == null) {
            shard.signalDone();
            bound = new BoundedReadIterator(iter, 0);
        }
        return bound;
    }

    private SamFileHeaderMerger CreateHeader() {
        // TODO: make extremely less horrible
        List<SAMFileReader> lst = GetReaderList();

        // now merge the headers
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(lst, SORT_ORDER);
        return headerMerger;
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
        BoundedReadIterator bound;// is this the first time we're doing this?
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

        // we're good, increment our read cout
        this.readsTaken += readCount;
        return new BoundedReadIterator(iter, readCount);

    }


    /**
     * unmapped reads.
     *
     * @param readCount how many reads to retrieve
     * @param iter      the iterator to use
     * @return the bounded iterator that you can use to get the intervaled reads from
     * @throws SimpleDataSourceLoadException
     */
    private BoundedReadIterator fastMappedReadSeek(long readCount, MergingSamRecordIterator2 iter) throws SimpleDataSourceLoadException {
        BoundedReadIterator bound;// is this the first time we're doing this?
        if (lastReadPos == null) {
            lastReadPos = new GenomeLoc(iter.getMergedHeader().getSequenceDictionary().getSequence(0).getSequenceIndex(), 0, 0);
            iter.queryContained(lastReadPos.getContig(), 1, -1);
            bound = new BoundedReadIterator(iter, readCount);
            this.readsTaken = readCount;
        }
        // we're not at the beginning, not at the end, so we move forward with our ghastly plan...
        else {

            iter.queryContained(lastReadPos.getContig(), (int) lastReadPos.getStop(), -1);

            // move the number of reads we read from the last pos
            while (iter.hasNext() && this.readsAtLastPos > 0) {
                iter.next();
                --readsAtLastPos;
            }
            if (readsAtLastPos != 0) {
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
                        SamFileHeaderMerger mg = CreateHeader();
                        iter = new MergingSamRecordIterator2(mg);
                        iter.queryContained(lastReadPos.getContig(), 1, Integer.MAX_VALUE);
                        return new BoundedReadIterator(iter,readCount);
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
                throw new StingException("Danger");
                //return null;
            }
            bound = new BoundedReadIterator(iter, readCount);
        }


        // return the iterator
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
