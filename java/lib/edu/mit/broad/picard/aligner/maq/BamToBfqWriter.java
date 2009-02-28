/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.aligner.maq;

import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.filter.*;
import edu.mit.broad.picard.util.PeekableIterator;
import edu.mit.broad.picard.util.Log;
import edu.mit.broad.picard.sam.ReservedTagConstants;

import java.io.File;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Arrays;

/**
 * Class to take unmapped reads in BAM file format and create Maq binary fastq format file(s) --
 * one or two of them, depending on whether it's a paired-end read.  This relies on the unmapped
 * BAM file having all paired reads together in order.
 */
public class BamToBfqWriter {

    private final File bamFile;
    private final String outputPrefix;
    private boolean pairedReads = false;
    private int wrote = 0;
    private int increment = 1;
    private int chunk = 0;
    private BinaryCodec codec1;
    private BinaryCodec codec2;
    private final Log log = Log.getInstance(BamToBfqWriter.class);

    /**
     * Constructor
     *
     * @param bamFile        the BAM file to read from
     * @param outputPrefix   the directory and file prefix for the binary fastq files
     * @param total          the total number of records that should be written, drawn evenly
     *                       from throughout the file (null for all).
     * @param chunk          the maximum number of records taht should be written to any one file
     * @param pairedReads    whether these reads are from  a paired-end run
     */
    public BamToBfqWriter(File bamFile, String outputPrefix, Integer total, Integer chunk, boolean pairedReads) {
        this.bamFile = bamFile;
        this.outputPrefix = outputPrefix;
        this.pairedReads = pairedReads;
        if (total != null) {
            double writeable = (double)countWritableRecords();
            this.increment = (int)Math.floor(writeable/total.doubleValue());
        }
        if (chunk != null) {
            this.chunk = chunk;
        }
    }

    /**
     * Constructor
     *
     * @param bamFile   the BAM file to read from
     * @param outputPrefix   the directory and file prefix for the binary fastq files
     * @param pairedReads    whether these reads are from  a paired-end run
     */
    public BamToBfqWriter(File bamFile, String outputPrefix, boolean pairedReads) {
        this(bamFile, outputPrefix, null, null, pairedReads);
    }

    /**
     * Writes the binary fastq file(s) to the output directory
     */
    public void writeBfqFiles() {

        Iterator<SAMRecord> iterator = (new SAMFileReader(IoUtil.openFileForReading(this.bamFile))).iterator();

        // Filter out noise reads and reads that fail the quality filter
        TagFilter tagFilter = new TagFilter(ReservedTagConstants.XN, 1);
        FailsVendorReadQualityFilter qualityFilter = new FailsVendorReadQualityFilter();

        if (!pairedReads) {
            writeSingleEndBfqs(iterator, Arrays.asList(tagFilter, qualityFilter));
            codec1.close();
        }
        else {
            writePairedEndBfqs(iterator, tagFilter, qualityFilter);
            codec1.close();
            codec2.close();
        }
        log.info("Wrote " + wrote + " bfq records.");

    }

    /**
     * Path for writing bfqs for paired-end reads
     *
     * @param iterator      the iterator witht he SAM Records to write
     * @param tagFilter     the filter for noise reads
     * @param qualityFilter the filter for PF reads
     */
    private void writePairedEndBfqs(Iterator<SAMRecord> iterator, TagFilter tagFilter,
                                    FailsVendorReadQualityFilter qualityFilter) {
        // Open the codecs for writing
        int fileIndex = 0;
        initializeNextBfqFiles(fileIndex++);

        int records = 0;

        while (iterator.hasNext()) {
            SAMRecord first = iterator.next();
            if (!iterator.hasNext()) {
                throw new PicardException("Mismatched number of records in " + this.bamFile.getAbsolutePath());
            }
            SAMRecord second = iterator.next();
            if (!second.getReadName().equals(first.getReadName()) ||
                first.getFirstOfPairFlag() == second.getFirstOfPairFlag()) {
                throw new PicardException("Unmatched read pairs in " + this.bamFile.getAbsolutePath() +
                    ": " + first.getReadName() + ", " + second.getReadName() + ".");
            }

            // If both are noise reads, filter them out
            if (tagFilter.filterOut(first) && tagFilter.filterOut(second))  {
                // skip it
            }
            // If either fails to pass filter, then exclude them as well
            else if (qualityFilter.filterOut(first) || qualityFilter.filterOut(second)) {
                // skip it
            }
            // Otherwise, write them out
            else {
                records++;
                if (records % increment == 0) {
                    first.setReadName(first.getReadName() + "#0/1");
                    writeFastqRecord(first.getFirstOfPairFlag() ? codec1 : codec2, first);
                    second.setReadName(second.getReadName() + "#0/2");
                    writeFastqRecord(second.getFirstOfPairFlag() ? codec1 : codec2, second);
                    wrote++;
                    if (wrote % 1000000 == 0) {
                        log.info(wrote + " records written.");
                    }
                    if (chunk > 0 && wrote % chunk == 0) {
                        initializeNextBfqFiles(fileIndex++);
                    }
                }
            }
        }
    }

    /**
     * Path for writing bfqs for single-end reads
     *
     * @param iterator  the iterator witht he SAM Records to write
     * @param filters   the list of filters to be applied
     */
    private void writeSingleEndBfqs(Iterator<SAMRecord> iterator, List<SamRecordFilter> filters) {

        // Open the codecs for writing
        int fileIndex = 0;
        initializeNextBfqFiles(fileIndex++);

        int records = 0;

        FilteringIterator it = new FilteringIterator(iterator, new AggregateFilter(filters));
        while (it.hasNext()) {
            SAMRecord record = it.next();
            records++;
            if (records % increment == 0) {

                writeFastqRecord(codec1, record);
                wrote++;
                if (wrote % 1000000 == 0) {
                    log.info(wrote + " records processed.");
                }
                if (chunk > 0 && wrote % chunk == 0) {
                    initializeNextBfqFiles(fileIndex++);
                }
            }
        }
    }

    /**
     * Closes any the open bfq file(s), if any, and opens the new one(s)
     *
     * @param fileIndex the index (counter) of the files to write
     */
    private void initializeNextBfqFiles(int fileIndex) {
        // Close the codecs if they were writing before
        if (codec1 != null) {
            codec1.close();
            if (pairedReads) {
                codec2.close();
            }
        }

        // Open new file, using the fileIndex.
        File bfq1 = getOutputFile(this.outputPrefix , 1, fileIndex);
        codec1 = new BinaryCodec(IoUtil.openFileForWriting(bfq1));
        log.info("Now writing to file " + bfq1.getAbsolutePath());
        if (pairedReads) {
            File bfq2 = getOutputFile(this.outputPrefix , 2, fileIndex);
            codec2 = new BinaryCodec(IoUtil.openFileForWriting(bfq2));
            log.info("Now writing to file " + bfq2.getAbsolutePath());
        }
    }

    /**
     * Writes out a SAMRecord in Maq fastq format
     *
     * @param codec the code to write to
     * @param rec   the SAMRecord to write
     */
    private void writeFastqRecord(BinaryCodec codec, SAMRecord rec) {

        // Writes the length of the read name and then the name (null-terminated)
        codec.writeString(rec.getReadName(), true, true);

        char seqs[] = rec.getReadString().toCharArray();
        char quals[] = rec.getBaseQualityString().toCharArray();
        
        // Write the length of the sequence
        codec.writeInt(seqs.length);

        // Calculate and write the sequence and qualities
        byte seqsAndQuals[] = new byte[seqs.length];

        for (int i = 0; i < seqs.length; i++) {
            int quality = Math.min(quals[i]-33, 63);
            int base;
            switch(seqs[i]) {
                case 'A':
                case 'a':
                    base = 0;
                    break;
                case 'C':
                case 'c':
                    base = 1;
                    break;
                case 'G':
                case 'g':
                    base = 2;
                    break;
                case 'T':
                case 't':
                    base = 3;
                    break;
                case 'N':
                case 'n':
                case '.':
                    base = 0;
                    quality = 0;
                    break;
                default:
                    throw new PicardException("Unknown base when writing bfq file: " + seqs[i]);
            }
            seqsAndQuals[i] = (byte) (base << 6 | quality);
        }
        codec.writeBytes(seqsAndQuals);
    }

    private int countWritableRecords() {
        int count = 0;
        PeekableIterator<SAMRecord> it = new PeekableIterator<SAMRecord>((new SAMFileReader(IoUtil.openFileForReading(this.bamFile))).iterator());
        if (!this.pairedReads) {
            // Filter out noise reads and reads that fail the quality filter
            List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
            filters.add(new TagFilter(ReservedTagConstants.XN, 1));
            filters.add(new FailsVendorReadQualityFilter());
            FilteringIterator itr = new FilteringIterator(it, new AggregateFilter(filters));
            while (itr.hasNext()) {
                itr.next();
                count++;
            }
        }
        else {
            while (it.hasNext()) {
                SAMRecord first = it.next();
                SAMRecord second = it.next();
                // If both are noise reads, filter them out
                if (first.getAttribute(ReservedTagConstants.XN) != null &&
                    second.getAttribute(ReservedTagConstants.XN) != null)  {
                    // skip it
                }
                // If either fails to pass filter, then exclude them as well
                else if (first.getReadFailsVendorQualityCheckFlag() || second.getReadFailsVendorQualityCheckFlag() ) {
                    // skip it
                }
                // Otherwise, write them out
                else {
                    count++;
                }
            }
        }
        it.close();
        return count;
    }

    /**
     * Constructs the name for the output file and returns the file
     *
     * @param outputPrefix        the directory and file prefix for the output bfq file
     * @param read                whether this is the file for the first or second read
     * @return                    a new File object for the bfq file.
     */
    private File getOutputFile(String outputPrefix, int read, int index) {
        File result = new File(outputPrefix + "." + index + "." + read + ".bfq");
        IoUtil.assertFileIsWritable(result);
        return result;
    }

}
