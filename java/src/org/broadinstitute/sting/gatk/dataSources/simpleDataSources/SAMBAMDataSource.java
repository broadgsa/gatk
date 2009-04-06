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

    // our sam file readers
    private final ArrayList<SAMFileReader> readers = new ArrayList<SAMFileReader>();

    // do we care that the SAM files respect the sort order.
    private boolean matchedSortOrders = true;

    // our merged sam iterator for spliting up the files
    MergingSamRecordIterator2 mergeIterator;

    // are we set to locus mode or read mode for dividing
    private boolean locusMode = true;

    // How strict should we be with SAM/BAM parsing?
    protected SAMFileReader.ValidationStringency strictness = SAMFileReader.ValidationStringency.STRICT;

    /**
     * constructor, given a single sam file
     *
     * @param samFiles the list of sam files
     */
    public SAMBAMDataSource(List<String> samFiles) throws SimpleDataSourceLoadException {
        List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        for (String fileName : samFiles) {
            File smFile = new File(fileName);
            if (!smFile.canRead()) {
                throw new SimpleDataSourceLoadException("SAMBAMDataSource: Unable to load file: " + fileName);
            }
            SAMFileReader reader = initializeSAMFile(smFile);
            if (reader != null) {
                readers.add(reader);
            }
        }

        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SORT_ORDER);
        this.mergeIterator = new MergingSamRecordIterator2(headerMerger);
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
     * getQueryRegionIterator
     * </p>
     *
     * @param location the genome location to extract data for
     * @return an iterator for that region
     */
    public MergingSamRecordIterator2 seek(GenomeLoc location) {
        MergingSamRecordIterator2 iter = new MergingSamRecordIterator2(this.mergeIterator);
        if (locusMode) {
            iter.query(location.getContig(), (int) location.getStart(), (int) location.getStop(), true);
        } else {
            iter.queryContained(location.getContig(), (int) location.getStart(), (int) location.getStop());
        }
        return iter;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
