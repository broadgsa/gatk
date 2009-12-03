package org.broadinstitute.sting.gatk.walkers;

/*  Motivation: BAM files created by samtools which are sorted often don't
 *   have the sort order set to 'coordinate' in the BAM header (instead it's
 *   marked as 'unsorted'.  Because BAMs are binary files, there's no way to
 *   easily change the tag.
 *
 *  This walker rewrites a BAM file so that the output is identical to the input
 *   except that the sort order tag is set to 'coordinate'.
 *
 *  Important note: to run properly in the GATK, the '-U' flag must be used so that
 *   the input BAM file not be rejected by the system.
 */

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriterFactory;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Oct 9, 2009
 * Time: 2:21:08 PM
 */
public class FixBAMSortOrderTag extends ReadWalker<SAMRecord, SAMFileWriter> {

    @Argument(required=true, shortName="out_bam", doc="The samfile to output to")
    public File SAM_FILE_OUTPUT_LOCATION;

    @Argument(required=false, shortName="sortorder", doc="the sort order to emit in")
    public SAMFileHeader.SortOrder SORT_ORDER=SAMFileHeader.SortOrder.coordinate;

    @Override
    public SAMRecord map(char[] ref, SAMRecord read) {
        return read;
    }

    @Override
    public SAMFileWriter reduceInit() {
        SAMFileHeader header = GenomeAnalysisEngine.instance.getSAMFileHeader();
        header.setSortOrder(SORT_ORDER);
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        return factory.makeBAMWriter(header,false,SAM_FILE_OUTPUT_LOCATION);
    }

    @Override
    public SAMFileWriter reduce(SAMRecord value, SAMFileWriter sum) {
        sum.addAlignment(value);
        return sum;
    }

    public void onTraversalDone(SAMFileWriter result) {
        result.close();
    }
}
