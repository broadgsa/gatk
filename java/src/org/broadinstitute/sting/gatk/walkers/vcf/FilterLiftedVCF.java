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

import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broad.tribble.vcf.VCFRecord;
import org.broad.tribble.vcf.VCFCodec;

import java.util.List;

/**
 * Filters a lifted-over VCF file for ref bases that have been changed.
 */
@Requires(value={},referenceMetaData=@RMD(name="vcf",type= VCFCodec.class))
public class FilterLiftedVCF extends RodWalker<Integer, Integer> {

    private VCFWriter writer;

    private long failedLocs = 0, totalLocs = 0;

    public void initialize() {}

    private void filterAndWrite(char ref, VCFRecord record) {

        totalLocs++;

        char recordRef = record.getReference().charAt(0);

        if ( recordRef != ref ) {

            // is it reverse complemented?
            if ( BaseUtils.simpleComplement(recordRef) == ref ) {
                record.setFilterString(record.isFiltered() ? String.format("%s;LiftoverToReverseComplementRefBase", record.getFilterString()) : "LiftoverToNewBase");
                failedLocs++;
            } else {
                record.setFilterString(record.isFiltered() ? String.format("%s;LiftoverToDifferentRefBase", record.getFilterString()) : "LiftoverToNewBase");
                failedLocs++;
            }
        }

        if ( writer == null ) {
            writer = new VCFWriter(out);
            writer.writeHeader(record.getHeader());
        }
        writer.addRecord(record);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        List<Object> rods = tracker.getReferenceMetaData("vcf");

        for ( Object rod : rods )
            filterAndWrite(ref.getBase(), (VCFRecord)rod);

        return 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return 0; }

    public void onTraversalDone(Integer result) {
        if ( writer != null )
            writer.close();
        System.out.println("Filtered " + failedLocs + " records out of " + totalLocs + " total records.");
    }
}