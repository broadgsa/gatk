package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
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
public class UpdateSAMSortOrderWalker extends ReadWalker<SAMRecord, SAMFileWriter> {

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
