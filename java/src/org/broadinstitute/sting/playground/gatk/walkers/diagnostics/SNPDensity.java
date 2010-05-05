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

package org.broadinstitute.sting.playground.gatk.walkers.diagnostics;

import org.broad.tribble.vcf.VCFCodec;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.commandline.Argument;

/**
 * Computes the density of SNPs passing and failing filters in intervals on the genome and emits a table for display
 */
@By(DataSource.REFERENCE)
@Requires(value={},referenceMetaData=@RMD(name="eval",type= VCFCodec.class))
public class SNPDensity extends RefWalker<Pair<VariantContext, GenomeLoc>, SNPDensity.Counter> {
    @Argument(fullName="granularity", shortName="granularity", doc="", required=false)
    private int granularity = 1000000;

    public void initialize() {
        out.printf("chr middlePos   linearPos   nSNPs   nSNPsFiltered   unfiltered.density filtered.density%n");
    }

    public class Counter {
        GenomeLoc firstLoc = null;
        long linearOffset = 0;
        int nSNPsCalled = 0;
        int nSNPsFiltered = 0;

        public Counter(Long linearOffset) {
            this.linearOffset = linearOffset;
            //System.out.printf("linear offset %d%n", linearOffset);
        }
    }

    public Pair<VariantContext, GenomeLoc> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContext vc = null;

        VCFRecord vcf = tracker.lookup("eval",VCFRecord.class);
        if (vcf != null)
            vc = VariantContextAdaptors.toVariantContext("eval", vcf);
        return new Pair<VariantContext, GenomeLoc>(vc, context.getLocation());
    }

    public Counter reduceInit() {
        return new Counter(0L);
    }

    private void printLine(Counter sum) {
        long offset = granularity / 2 - 1;
        long chrOffset = sum.firstLoc.getStart() + offset;
        out.printf("%s  %d  %d    %d    %d  %.2e    %.2e%n",
                sum.firstLoc.getContig(),
                chrOffset,
                sum.linearOffset + offset,
                sum.nSNPsCalled, sum.nSNPsFiltered,
                (1.0 * sum.nSNPsCalled) / granularity, (1.0 * sum.nSNPsFiltered) / granularity);
    }

    public Counter reduce(Pair<VariantContext, GenomeLoc> p, Counter sum) {
        if ( p == null )
            return sum;

//        System.out.printf("%s %s %d%n", c.getLocation(), sum.firstLoc, sum.nSNPsSeen);
        VariantContext c = p.getFirst();
        GenomeLoc loc = p.getSecond();

        if ( sum.firstLoc != null ) {
            long dist = loc.distance(sum.firstLoc);
//            System.out.printf("  dist = %d%n", dist);
            if ( dist > granularity ) {
                printLine(sum);
                sum = new Counter(sum.linearOffset + granularity);
            }
        }

        if ( sum.firstLoc == null ) sum.firstLoc = loc;

        sum.nSNPsCalled         += c != null && c.isNotFiltered() ? 1 : 0;
        sum.nSNPsFiltered += c != null && c.isFiltered() ? 1 : 0;

        return sum;
    }

    public void onTraversalDone(Counter sum) {
        printLine(sum);
    }
}