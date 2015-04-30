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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Filters a lifted-over VCF file for reference bases that have been changed
 *
 * <p>"Lifting over" variants means adjusting variant calls from one reference to another. Specifically, the process
 * adjusts the position of the call to match the corresponding position on the target reference. For example, if you
 * have variants called from reads aligned to the hg19 reference, and you want to compare them to calls made based on
 * the b37 reference, you need to liftover one of the callsets to the other reference.</p>
 *
 * <p>This tool is intended to be the second of two processing steps for the liftover process. The first step is to
 * run LiftoverVariants on your VCF file. The second step is to run FilterLiftedVariants on the output of
 * LiftoverVariants. This will produce valid well-behaved VCF files, where you'll see that the contig names in the
 * header have all been correctly replaced.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A lifted-over variant call set to filter.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * The filtered call set.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T FilterLiftedVariants \
 *   -R reference.fasta \
 *   -V liftedover_input.vcf \
 *   -o filtered_output.vcf
 * </pre>
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

    /**
     * Determines whether records should be filtered; if not, writes them to the output
     *
     * @param ref   the reference context
     * @param vc    the VariantContext to process
     * @return true if the record is not filtered, false otherwise
     */
    protected boolean filterOrWrite(final byte[] ref, final VariantContext vc) {
	if ( ref == null ) throw new IllegalArgumentException("Cannot filter based on a null reference array");
	if ( vc == null ) throw new IllegalArgumentException("Cannot filter a null Variant Context");

        totalLocs++;

        boolean filter = false;
        final byte[] recordRef = vc.getReference().getBases();

        // this can happen for records that get placed at the ends of chromosomes
        if ( recordRef.length > ref.length ) {
            filter = true;
        } else {
            for (int i = 0; i < recordRef.length && i < MAX_VARIANT_SIZE; i++) {
                if ( recordRef[i] != ref[i] ) {
                    filter = true;
                    break;
                }
            }
        }

        if ( filter )
            failedLocs++;
        else
            writer.add(vc);

        return !filter;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        final Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
        for ( final VariantContext vc : VCs )
            filterOrWrite(ref.getBases(), vc);

        return 0;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return 0; }

    public void onTraversalDone(Integer result) {
        System.out.println("Filtered " + failedLocs + " records out of " + totalLocs + " total records.");
    }
}
