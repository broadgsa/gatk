/*
* Copyright 2012-2015 Broad Institute, Inc.
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

import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.util.*;

/**
 * Simplify multi-nucleotide variants (MNPs) into more basic/primitive alleles.
 *
 * <p>This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component
 * part (A-T and A->G).</p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant set with any type of alleles.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A VCF with alleles broken into primitive types.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T VariantsToAllelicPrimitives \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>For now this tool modifies only multi-nucleotide polymorphisms (MNPs) and leaves SNPs, indels, and
 * complex substitutions as is, although one day it may be extended to handle the complex substitution case.</li>
 *     <li>This tool modifies only bi-allelic variants.</li>
 * </ul>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class VariantsToAllelicPrimitives extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter baseWriter = null;

    private VariantContextWriter vcfWriter;

    public void initialize() {
        final String trackName = variantCollection.variants.getName();
        final Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(trackName));

        final Map<String, VCFHeader> vcfHeaders = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));
        final Set<VCFHeaderLine> headerLines = vcfHeaders.get(trackName).getMetaDataInSortedOrder();

        baseWriter.writeHeader(new VCFHeader(headerLines, samples));

        vcfWriter = VariantContextWriterFactory.sortOnTheFly(baseWriter, 200);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        final Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());

        int changedSites = 0;
        for ( final VariantContext vc : VCs )
            changedSites += writeVariants(vc);

        return changedSites;
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        System.out.println(result + " MNPs were broken up into primitives");
        vcfWriter.close();
    }

    @Requires("vc != null")
    private int writeVariants(final VariantContext vc) {
        // for now, we modify only bi-allelic MNPs; update docs above if this changes
        if ( vc.isBiallelic() && vc.isMNP() ) {
            for ( final VariantContext splitVC : GATKVariantContextUtils.splitIntoPrimitiveAlleles(vc) )
                vcfWriter.add(splitVC);
            return 1;
        } else {
            vcfWriter.add(vc);
            return 0;
        }
    }
}
