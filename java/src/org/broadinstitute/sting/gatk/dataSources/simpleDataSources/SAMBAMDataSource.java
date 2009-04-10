package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import edu.mit.broad.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.utils.GenomeLoc;

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
public class SAMBAMDataSource implements SimpleDataSource {
    /** our SAM data files */
    private final SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMBAMDataSource.class);

    // are we set to locus mode or read mode for dividing
    private boolean locusMode = true;

    // How strict should we be with SAM/BAM parsing?
    protected SAMFileReader.ValidationStringency strictness = SAMFileReader.ValidationStringency.STRICT;

    // our list of readers
    private final List<File> samFileList = new ArrayList<File>();

    /**
     * constructor, given a single sam file
     *
     * @param samFiles the list of sam files
     */
    public SAMBAMDataSource(List<String> samFiles) throws SimpleDataSourceLoadException {
        for (String fileName : samFiles) {
            File smFile = new File(fileName);
            if (!smFile.canRead()) {
                throw new SimpleDataSourceLoadException("SAMBAMDataSource: Unable to load file: " + fileName);
            }
            samFileList.add(smFile);

        }

        //SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(samFileList, SORT_ORDER);
    }


    protected SAMFileReader initializeSAMFile(final File samFile) {
        if (samFile.toString().endsWith(".list")) {
            return null;
        } else {
            SAMFileReader samReader = new SAMFileReader(samFile, true);
            samReader.setValidationStringency(strictness);

            final SAMFileHeader header = samReader.getFileHeader();
            logger.info(String.format("Sort order is: " + header.getSortOrder()));

            return samReader;
        }
    }


    /**
     * set the mode to by loci, which let's you duplicate reads, but never at a single
     * locus, or false for read mode where no read is seen twice.
     *
     * @param tr true if by loci, false if by read
     */
    public void setToByLociMode(boolean tr) {
        locusMode = tr;
    }

    /**
     * <p>
     * seek
     * </p>
     *
     * @param location the genome location to extract data for
     * @return an iterator for that region
     */
    public MergingSamRecordIterator2 seek(GenomeLoc location) throws SimpleDataSourceLoadException {

        // right now this is pretty damn heavy, it copies the file list into a reader list every time
        List<SAMFileReader> lst = new ArrayList<SAMFileReader>();
        for (File f : this.samFileList) {
            SAMFileReader reader = initializeSAMFile(f);
            if (reader == null) {
                throw new SimpleDataSourceLoadException("SAMBAMDataSource: Unable to load file: " + f);
            }
            lst.add(reader);
        }

        // now merge the headers
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(lst, SORT_ORDER);

        // make a merging iterator for this record
        MergingSamRecordIterator2 iter =  new MergingSamRecordIterator2(headerMerger);


        // we do different things for locus and read modes
        if (locusMode) {
            iter.queryOverlapping(location.getContig(), (int) location.getStart(), (int) location.getStop());
        } else {
            iter.queryContained(location.getContig(), (int) location.getStart(), (int) location.getStop());
        }

        // return the iterator
        return iter;
    }

    

}
