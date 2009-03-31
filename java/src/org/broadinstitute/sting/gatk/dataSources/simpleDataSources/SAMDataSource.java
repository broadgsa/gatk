package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.sam.MergingSamRecordIterator;
import edu.mit.broad.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;

import java.io.File;
import java.io.FileNotFoundException;
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

    // our sam file readers
    private final ArrayList<SAMFileReader> readers = new ArrayList<SAMFileReader>();

    // do we care that the SAM files respect the sort order.
    private boolean matchedSortOrders = true;

    // our record iterator, we use it to iterate over all the reads
    private MergingSamRecordIterator iterator = null;

    // we may want to write out the file
    private SAMFileWriter out = null;

    // are we set to locus mode or read mode for dividing
    private boolean locusMode = true;

    /**
     * constructor for multiple sam files
     *
     * @param samfiles
     */
    public SAMDataSource(ArrayList<String> samfiles) throws FileNotFoundException {
        loadFiles(samfiles);
    }

    private void loadFiles(ArrayList<String> samfiles) throws FileNotFoundException {
        // verify the list passed to the class
        ArrayList<File> INPUT = new ArrayList<File>();
        for (String check : samfiles) {
            File nf = new File(check);
            if (!nf.exists()) {
                throw new FileNotFoundException(check + " doesn't exist");
            }
        }


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
        if (matchedSortOrders || SORT_ORDER == SAMFileHeader.SortOrder.unsorted) {
            SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SORT_ORDER);
            iterator = new MergingSamRecordIterator(headerMerger);

        } else {
            SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SAMFileHeader.SortOrder.unsorted);
            iterator = new MergingSamRecordIterator(headerMerger);
            SAMFileHeader header = headerMerger.getMergedHeader();
            header.setSortOrder(SORT_ORDER);

        }
    }

    /**
     * constructor, given a single sam file
     *
     * @param samFile
     */
    public SAMDataSource(String samFile) throws FileNotFoundException {
        ArrayList<String> samfiles = new ArrayList<String>();
        samfiles.add(samFile);
        loadFiles(samfiles);
    }

    /**
     * Chunk the sam file at appropriate locations, given the chunk count
     *
     * @param chunkCount
     * @return
     */
    public void chunk(int chunkCount) {

    }

    /** set this source to divide on reads */
    public void setToReadMode() {
        locusMode = true;
    }
}
