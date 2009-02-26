/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

import edu.mit.broad.sam.*;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.util.Log;
import edu.mit.broad.picard.filter.AggregateFilter;
import edu.mit.broad.picard.filter.SamRecordFilter;
import edu.mit.broad.picard.filter.SolexaNoiseFilter;
import edu.mit.broad.picard.sam.ReservedTagConstants;

import java.io.File;
import java.util.*;

/**
 * Writes the data from a BustardFileParser to a BAM file
 */
public class BustardToSamWriter {

    private final BustardFileParser parser;
    private SAMFileWriter writer;
    private final File outputFile;
    private AggregateFilter filters;
    private int recordsWritten = 0;
    private Log log = Log.getInstance(BustardToSamWriter.class);

    /**
     * Constructor
     *
     * @param parser            The parser for the Bustard data
     * @param outputDirectory   The directory in which to write the BAM file
     * @param flowcell          The flowcell from which the data is drawn
     */
    public BustardToSamWriter(BustardFileParser parser, File outputDirectory, String flowcell) {
        this.parser = parser;
        this.outputFile = getOutputFile(outputDirectory, flowcell);
        initializeFilters();
    }

    /**
     * Alternate constructor for testing
     *
     * @param parser            The parser for the Bustard data
     * @param outputFile   The directory in which to write the BAM file
     */
    BustardToSamWriter(BustardFileParser parser, File outputFile) {
        this.parser = parser;
        this.outputFile = outputFile;
        initializeFilters();
    }

    private void initializeFilters() {
        filters = new AggregateFilter(Arrays.asList(
            (SamRecordFilter)new SolexaNoiseFilter()
        ));
    }


    /**
     * Writes all data from the BustardFileParser to a BAM file
     */
    public void writeBamFile() {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        writer = new SAMFileWriterFactory().makeBAMWriter(header, false, outputFile);

        while (parser.hasNext()) {
            BustardReadData brd = parser.next();

            SAMRecord sam = createSamRecord(brd, true);
            writer.addAlignment(sam);
            this.recordsWritten++;

            if (parser.isPairedEnd()) {
                SAMRecord sam2 = createSamRecord(brd, false);
                writer.addAlignment(sam2);
                this.recordsWritten++;
            }

        }
        writer.close();

        log.info("Wrote " + this.recordsWritten + " read records to BAM file " +
                this.outputFile.getAbsolutePath());
    }

    /**
     * Creates a SAMRecord from Bustard data
     *
     * @param brd           The BustardReadData to use in populating the SAMRecord
     * @param isFirstRead   whether this is the first read of a pair
     * @return SAMRecord    fully populated SAMRecord
     */
    private SAMRecord createSamRecord(BustardReadData brd, boolean isFirstRead) {
        SAMRecord sam = new SAMRecord();
        sam.setReadName(brd.getReadName());
        sam.setReadString(isFirstRead ? brd.getFirstReadSequence() : brd.getSecondReadSequence());
        sam.setBaseQualityString(isFirstRead ? brd.getFirstReadPhredQualities() : brd.getSecondReadPhredQualities());

        // Flag values
        sam.setReadPairedFlag(parser.isPairedEnd());
        sam.setReadUmappedFlag(true);
        sam.setReadFailsVendorQualityCheckFlag(!brd.isPf());
        sam.setMateUnmappedFlag(true);
        if (parser.isPairedEnd()) {
            sam.setFirstOfPairFlag(isFirstRead);
            sam.setSecondOfPairFlag(!isFirstRead);
        }

        if (filters.filterOut(sam)) {
            sam.setAttribute(ReservedTagConstants.XN, 1);
        }
        return sam;
    }

    /**
     * Constructs the name for the output file, determines whether it is writeable,
     * and returns the file
     *
     * @param outputDirectory     the directory in which to write the BAM file
     * @param flowcell            the flowcell from which the data is drawn
     * @return                    a new File object for the BAM file.
     */
    private File getOutputFile(File outputDirectory, String flowcell) {
        File result = new File(outputDirectory.getAbsolutePath() + "/" +
                flowcell + "." + parser.getLane() + ".unmapped.bam");
        IoUtil.assertFileIsWritable(result);
        return result;
    }
}
