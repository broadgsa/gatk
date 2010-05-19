/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

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
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.commandline.Argument;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriterFactory;

import java.io.File;

/**
 * Fixes slightly corrupted BAM files by rewriting the input BAM file, altering
 * the header by changing the sort order tag (SO) to coordinate sort order.  Will NOT
 * verify the contents of the file to ensure that the data is actually in coordinate sorted
 * order. 
 */
public class FixBAMSortOrderTag extends ReadWalker<SAMRecord, SAMFileWriter> {

    @Argument(required=true, shortName="out_bam", doc="The samfile to output to")
    public File SAM_FILE_OUTPUT_LOCATION;

    @Argument(required=false, shortName="sortorder", doc="the sort order to emit in")
    public SAMFileHeader.SortOrder SORT_ORDER=SAMFileHeader.SortOrder.coordinate;

    @Override
    public SAMRecord map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
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
