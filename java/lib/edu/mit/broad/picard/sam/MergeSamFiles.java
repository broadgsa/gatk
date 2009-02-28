/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever.
* Neither the Broad Institute nor MIT can be responsible for its use, misuse, or
* functionality.
*/
package edu.mit.broad.picard.sam;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.sam.SAMFileHeader;
import static edu.mit.broad.sam.SAMFileHeader.SortOrder;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMFileWriter;
import edu.mit.broad.sam.SAMFileWriterFactory;
import edu.mit.broad.sam.SAMRecord;

/**
 * Reads a SAM or BAM file and combines the output to one file
 *
 * @author Dave Tefft
 */
public class MergeSamFiles extends CommandLineProgram {
    // Usage and parameters
    @Usage(programVersion="1.0")
    public String USAGE = "Merges multiple SAM/BAM files into one file.\n";

    @Option(shortName="I", doc="SAM or BAM input file", minElements=1)
    public List<File> INPUT = new ArrayList<File>();

    @Option(shortName="O", doc="SAM or BAM file to write merged result to")
    public File OUTPUT;

    @Option(shortName="SO", doc="Sort order of output file", optional=true)
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new MergeSamFiles().instanceMain(argv));
    }

    /** Combines multiple SAM/BAM files into one. */
    @Override
	protected int doWork() {
        boolean matchedSortOrders = true;

        // Open the files for reading and writing
        List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        for (File inFile : INPUT) {
            IoUtil.assertFileIsReadable(inFile);
            SAMFileReader in = new SAMFileReader(inFile);
            readers.add(in);
            matchedSortOrders = matchedSortOrders && in.getFileHeader().getSortOrder() == SORT_ORDER;
        }

        // If all the input sort orders match the output sort order then just merge them and
        // write on the fly, otherwise setup to merge and sort before writing out the final file
        IoUtil.assertFileIsWritable(OUTPUT);
        MergingSamRecordIterator iterator = null;
        SAMFileWriter out = null;

        if (matchedSortOrders || SORT_ORDER == SortOrder.unsorted) {
            SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SORT_ORDER);
            iterator = new MergingSamRecordIterator(headerMerger);
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerMerger.getMergedHeader(), true, OUTPUT);
        }
        else {
            SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SortOrder.unsorted);
            iterator = new MergingSamRecordIterator(headerMerger);
            SAMFileHeader header = headerMerger.getMergedHeader();
            header.setSortOrder(SORT_ORDER);
            out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);
        }

        // Lastly loop through and write out the records
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            out.addAlignment(record);
        }

        out.close();
        return 0;
    }

}