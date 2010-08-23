/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.StandardVCFWriter;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.commandline.Output;
import org.broad.tribble.vcf.VCFHeader;

import java.util.*;
import java.io.PrintStream;

/**
 * Filters a lifted-over VCF file for ref bases that have been changed.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant",type= ReferenceOrderedDatum.class))
public class FilterLiftedVariants extends RodWalker<Integer, Integer> {
    @Output
    PrintStream out;

    private VCFWriter writer;

    private long failedLocs = 0, totalLocs = 0;

    public void initialize() {
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList("variant"));
        Map<String, VCFHeader> vcfHeaders = VCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList("variant"));

        writer = new StandardVCFWriter(out);
        final VCFHeader vcfHeader = new VCFHeader(vcfHeaders.containsKey("variant") ? vcfHeaders.get("variant").getMetaData() : null, samples);
        writer.writeHeader(vcfHeader);
    }

    private void filterAndWrite(byte ref, VariantContext vc) {

        totalLocs++;

        byte recordRef = vc.getReference().getBases()[0];

        if ( recordRef != ref ) {
            failedLocs++;
        } else {
            writer.add(vc, ref);
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getVariantContexts(ref, "variant", null, context.getLocation(), true, false);
        for ( VariantContext vc : VCs )
            filterAndWrite(ref.getBase(), vc);

        return 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return 0; }

    public void onTraversalDone(Integer result) {
        writer.close();
        System.out.println("Filtered " + failedLocs + " records out of " + totalLocs + " total records.");
    }
}