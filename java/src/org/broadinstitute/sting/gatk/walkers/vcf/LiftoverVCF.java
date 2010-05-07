/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.vcf;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broad.tribble.vcf.VCFRecord;
import org.broad.tribble.vcf.VCFCodec;

import java.io.File;
import java.util.List;

import net.sf.picard.liftover.LiftOver;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;

/**
 * Lifts a VCF file over from one build to another.  Note that the resulting VCF could be mis-sorted.
 */
@Requires(value={},referenceMetaData=@RMD(name="vcf",type= VCFCodec.class))
public class LiftoverVCF extends RodWalker<Integer, Integer> {

    @Argument(fullName="chain", shortName="chain", doc="Chain file", required=true)
    protected File CHAIN = null;

    @Argument(fullName="newSequenceDictionary", shortName="dict", doc="Sequence .dict file for the new build", required=true)
    protected File NEW_SEQ_DICT = null;

    private VCFWriter writer;

    private LiftOver liftOver;

    private long successfulIntervals = 0, failedIntervals = 0;

    public void initialize() {
        liftOver = new LiftOver(CHAIN);
        liftOver.setLiftOverMinMatch(LiftOver.DEFAULT_LIFTOVER_MINMATCH);

        final SAMFileHeader toHeader = new SAMFileReader(NEW_SEQ_DICT).getFileHeader();
        liftOver.validateToSequences(toHeader.getSequenceDictionary());
    }

    private void convertAndWrite(VCFRecord record) {

        final Interval fromInterval = new Interval(record.getChr(), record.getStart(), record.getEnd());
        final Interval toInterval = liftOver.liftOver(fromInterval);

        if ( toInterval != null ) {
            record.setLocation(toInterval.getSequence(), toInterval.getStart());
            if ( writer == null ) {
                writer = new VCFWriter(out);
                writer.writeHeader(record.getHeader());
            }
            writer.addRecord(record);
            successfulIntervals++;
        } else {
            failedIntervals++;
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        List<Object> rods = tracker.getReferenceMetaData("vcf");

        for ( Object rod : rods )
            convertAndWrite((VCFRecord)rod);

        return 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return 0; }

    public void onTraversalDone(Integer result) {
        writer.close();
        System.out.println("Converted " + successfulIntervals + " intervals; failed to convert " + failedIntervals + " intervals.");
    }
}
