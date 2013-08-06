/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Filters a lifted-over VCF file for ref bases that have been changed.
 *
 * "Lifting over" variants means adjusting variant calls from one reference to another. Specifically, the process adjusts the position of the call to match the corresponding position on the target reference.
 * For example, if you have variants called from reads aligned to the hg19 reference, and you want to compare them to calls made based on the b37 reference, you need to liftover one of the callsets to the other reference.
 *
 * FilteredLiftedVariants is intended to be the second of two processing steps for the liftover process. The first step is to run LiftoverVariants on your VCF file.
 * The second step is to run FilterLiftedVariants on the output of LiftoverVariants. This will produce valid well-behaved VCF files, where you'll see that the contig names in the header have all been correctly replaced.
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=0,stop=100))
public class FilterLiftedVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    private static final int MAX_VARIANT_SIZE = 100;

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter writer = null;

    private long failedLocs = 0, totalLocs = 0;

    public void initialize() {
        String trackName = variantCollection.variants.getName();
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(trackName));
        Map<String, VCFHeader> vcfHeaders = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));

        final VCFHeader vcfHeader = new VCFHeader(vcfHeaders.containsKey(trackName) ? vcfHeaders.get(trackName).getMetaDataInSortedOrder() : Collections.<VCFHeaderLine>emptySet(), samples);
        writer.writeHeader(vcfHeader);
    }

    private void filterAndWrite(byte[] ref, VariantContext vc) {

        totalLocs++;

        boolean failed = false;
        byte[] recordRef = vc.getReference().getBases();
        for (int i = 0; i < recordRef.length && i < MAX_VARIANT_SIZE; i++) {
            if ( recordRef[i] != ref[i] ) {
                failed = true;
                break;
            }
        }

        if ( failed )
            failedLocs++;
        else
            writer.add(vc);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
        for ( VariantContext vc : VCs )
            filterAndWrite(ref.getBases(), vc);

        return 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return 0; }

    public void onTraversalDone(Integer result) {
        System.out.println("Filtered " + failedLocs + " records out of " + totalLocs + " total records.");
    }
}