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

import edu.mit.broad.sam.*;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.StringUtil;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.util.SamPairUtil;

import java.io.File;
import java.io.BufferedInputStream;
import java.util.*;

/**
 * Reads a Maq map file and returns an an iterator of SAMRecords and a populated header
 *
 * IMPORTANT!  Even though the reads in the map file are in coordinate order, this iterator
 * will not necessarily return them in that order.  For paired reads, both will be
 * returned only after *both* records have been seen.
 *
 * @author Kathleen Tibbetts
 */
public class MapFileIterator implements CloseableIterator<SAMRecord> {

    public static final int MATE_UNMAPPED_FLAG = 64;
    public static final int READ_UNMAPPED_FLAG = 192;

    private static final int READ_NAME_LENGTH = 36;
    private static final int MAP_FORMAT = -1;
    private static final int MAX_READ_LENGTH = 128;

    private static final byte ACGT[] = {'A', 'C', 'G', 'T'};

    public static final String PROGRAM_RECORD = "0";

    private long recordCount = 0L;
    private int recordsRead = 0;
    private BinaryCodec mapCodec;
    private final SAMFileHeader header;
    private final boolean pairedReads;
    private final boolean jumpingLibrary;
    private final List<SAMRecord> next = new ArrayList<SAMRecord>();
    private final Map<String, SAMRecord> pending = new HashMap<String, SAMRecord>();
    private final List<File> mapFiles = new LinkedList<File>();

    /**
     * Constructor.  Opens the map file, reads the record count and header from it,
     * creates the SAMFileHeader, and queues up the first read
     *
     * @param mapFile           The Maq map file to read
     * @param commandLine       The command line used to invoke Maq (for the header)
     * @param pairedReads       Whether this is a paired-end run
     */
    public MapFileIterator(String commandLine, boolean pairedReads, boolean jumpingLibrary, File... mapFile) {
        if (mapFile.length == 0) {
            throw new IllegalArgumentException("At least one map file must be provided.");
        }
        mapFiles.addAll(Arrays.asList(mapFile));

        this.pairedReads = pairedReads;
        this.jumpingLibrary = jumpingLibrary;

        header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMProgramRecord program = new SAMProgramRecord(PROGRAM_RECORD);
        program.setProgramVersion(MaqConstants.getProgramVersion());
        program.setCommandLine(commandLine);
        header.addProgramRecord(program);

        queueNextMapFile();
    }

    /**
     * Queues up the next map file
     *
     * @return  true if there's another map file to iterate over
     */
    private boolean queueNextMapFile() {

        // Close the old file
        if (mapCodec != null) {
            mapCodec.close();
        }

        // If there are no more map files, return fales
        if (mapFiles.size() == 0) {
            return false;
        }

        // Otherwise, open the next file and reset the recordsRead count
        mapCodec = new BinaryCodec(new BufferedInputStream(IoUtil.openFileForReading(mapFiles.remove(0))));
        int format = mapCodec.readInt();
        if (format != MAP_FORMAT) {
            mapCodec.close();
            throw new PicardException("Unrecognized Maq map file format: " + format);
        }
        recordsRead = 0;


        // Read the sequences out of the map file and set them on the header
        int sequenceCount = mapCodec.readInt();
        List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>();
        for (int i = 0; i < sequenceCount; i++) {
            int length = mapCodec.readInt();
            // Write the sequence name, trimming off the null terminator
            sequences.add(new SAMSequenceRecord(mapCodec.readString(length).substring(0, length-1)));
        }
        if (header.getSequences() == null || header.getSequences().size() == 0) {
            header.setSequences(sequences);
        }
        else {
            // TODO: Check that the sequences match and throw and exception if they don't
        }
        recordCount = mapCodec.readLong();

        readNext();
        return true;
    }

    /**
     * Closes the BinaryCodec reading the map file
     */
    public void close() {
        mapCodec.close();
    }

    /**
     * @return true if the iteration has more elements
     */
    public boolean hasNext() {
         return next.size() > 0;
    }

    /**
     * @return the next SAMRecord in the iteration
     * @throws NoSuchElementException if this is called when hasNext() returns false
     */
    public SAMRecord next() {
        if (!hasNext()) {
            throw new NoSuchElementException("No more elements in this iteration");
        }
        SAMRecord result = next.remove(0);
        readNext();
        return result;
    }

    /**
     * Reads the next element from the map file.  If we are done with it, we put it in the <code>next</code>
     * list; if we are waiting to see its mate, we put it in the <code>pending</code> map.  Calls itself
     * repeatedly until there is at least one element in <code>next</code>.
     */
    private void readNext() {

        // If there's already a record queued up, just return
        if (next.size() > 0) {
            return;
        }

        // If we've read all there is, then any remaining records in the pending map should be returned.  
        // If this is not a PE run, then the pending map will be empty and we're done.
        if (recordsRead == recordCount) {
            if (pending.size() > 0) {
                StringBuffer sb = new StringBuffer();
                for (String item : pending.keySet()) {
                    sb.append(item).append("\n");
                }
                throw new PicardException("MapFileIterator pending map should have been empty but contained " +
                        "the following records: " + sb.toString());
            }
            queueNextMapFile();
            return;
        }

        // Otherwise, we read until there is at least one record in the <code>next</code> list
        readMapRecord();
        if (next.size() == 0) {
            readNext();
        }
    }

    /**
     * Reads one record from the map file and throws it onto the pending map or the next list,
     * depending on whether we have already seen its mate
     */
    private void readMapRecord() {

        // Now that we've got all the data from the binary file, write a SAMRecord and add it to
        // the new BAM file
        SAMRecord record = new SAMRecord();
        record.setAttribute(SAMTag.PG.toString(), PROGRAM_RECORD);
        record.setReadPairedFlag(this.pairedReads);
        
        // the last base is the single-end mapping quality.
        byte seqsAndQuals[] = new byte[MAX_READ_LENGTH-1];
        mapCodec.readBytes(seqsAndQuals);

        byte singleEndMappingQualityOrIndelLength = mapCodec.readByte();

         // the length of the read
        int readLength = mapCodec.readUByte();
        setSeqsAndQuals(seqsAndQuals, readLength, record);

        // the final mapping quality (unless <code>flag</code> below is 130, then it is the
        //  position of the indel (or 0 if no indel)
        int mappingQuality = mapCodec.readUByte();

        // mismatches in the 28bp (higher 4 bits) and mismatches (lower 4 bits)
        mapCodec.readUByte();
        // sum of errors of the best hit
        mapCodec.readUByte();
        // counts of all 0- and 1-mismatch hits on the reference
        mapCodec.readUByte();
        mapCodec.readUByte();

        // A bitwise flag. See the Maq docs for its full meaning
        int flag = mapCodec.readUByte();

        // the lower mapQ of the two ends (equals map_qual if unpaired); if flag is 130: mapQ of its mate
        int altQual = mapCodec.readUByte();

        // Index of the sequence for this read
        record.setReferenceIndex((int)mapCodec.readUInt(), getHeader());

        // Start position and strand
        long pos = mapCodec.readUInt();
        int startPos = ((int)((pos>>1)& 0x7FFFFFFF)) + 1;
        record.setAlignmentStart(startPos);
        record.setReadNegativeStrandFlag((pos&1) == 1);

        // offset of the mate (zero if unpaired, or two ends mapped to different chr)
        mapCodec.readInt();

        // The read name
        byte nameBytes[] = new byte[READ_NAME_LENGTH];
        mapCodec.readBytes(nameBytes);
        String name = StringUtil.bytesToString(nameBytes).trim();
        if (this.pairedReads) {
            if (name.endsWith("/1")) {
                record.setFirstOfPairFlag(true);
                record.setSecondOfPairFlag(false);
            }
            else if (name.endsWith("/2")) {
                record.setFirstOfPairFlag(false);
                record.setSecondOfPairFlag(true);
            }
            else {
                throw new PicardException("Unrecognized ending for paired read name: " + name);                
            }
            name = name.substring(0, name.length()-2);
        }
        record.setReadName(name);


        if (flag != 130 || singleEndMappingQualityOrIndelLength == 0) { // No indel
            record.setCigarString(readLength + "M");
            record.setMappingQuality(mappingQuality);
        }
        else {  // Indel
            int indelPos = mappingQuality;
            String cigar = indelPos + "M" + Math.abs(singleEndMappingQualityOrIndelLength);
            int remaining = readLength - indelPos;
            if (singleEndMappingQualityOrIndelLength > 0) {
                cigar += "I" + (remaining - singleEndMappingQualityOrIndelLength) + "M";
            }
            else {
                cigar += "D" + remaining + "M";
            }
            record.setCigarString(cigar);
            // In the docs, it look like there is a mapping quality for the mate, do we use that?
            record.setMappingQuality(altQual);
        }

        if (!pairedReads) {
            record.setProperPairFlag(false);
            next.add(record);
        }
        else {
            record.setMateUnmappedFlag(flag == MATE_UNMAPPED_FLAG);
            SAMRecord mate = pending.remove(record.getReadName());

            if (mate != null) {
                boolean proper = SamPairUtil.isProperPair(record, mate, jumpingLibrary);
                record.setProperPairFlag(proper);
                mate.setProperPairFlag(proper);

                SamPairUtil.setMateInfo(record, mate);

                int insertSize = SamPairUtil.computeInsertSize(record, mate);
                record.setInferredInsertSize(insertSize);
                mate.setInferredInsertSize(insertSize);

                if (!mate.getMateUnmappedFlag()) {
                    next.add(record);
                }
                if (!record.getMateUnmappedFlag()) {
                    next.add(mate);
                }
            }
            else {
                pending.put(record.getReadName(), record);
            }
        }

        // TODO: Figure out what do do about noise reads long-term
        // Note that it is possible that we have lost a "Noise read" annotation at this point.  Since
        // we try to map a pair if only one of the reads is classified as "noise", then for any paired
        // reads where one was a noise read and one was not, we will lose the noise annotation on the
        // one noisy read.  We have discussed either re-doing the noise evaluation here, modifying the
        // read name to carry the noise flag through Maq, or changing what reads we give to Maq.

        recordsRead++;
        
    }

    /**
     * Decodes the sequence and the qualities and sets them on the SAMrecords
     *
     * @param seqsAndQuals  the list of seqs and quals
     * @param readLength    the length of the read
     * @param sam           the SAMRecord to populate
     */
    private void setSeqsAndQuals(byte seqsAndQuals[], int readLength, SAMRecord sam) {
        byte sequence[] = new byte[readLength];
        byte qualities[] = new byte[readLength];
        for (int i = 0; i < readLength; i++) {
            byte b = seqsAndQuals[i];
            qualities[i] = (byte)(b & 0x3F);
            if (b == 0) {
                sequence[i] = 'N';
            }
            else {
                sequence[i] = ACGT[(seqsAndQuals[i] >> 6) & 3];
            }
        }
        sam.setReadBases(sequence);
        sam.setBaseQualities(qualities);
    }

    /**
     * @throws UnsupportedOperationException -- not implemented
     */
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported in MapFileIterator");
    }

    public SAMFileHeader getHeader() { return header;  }
}
