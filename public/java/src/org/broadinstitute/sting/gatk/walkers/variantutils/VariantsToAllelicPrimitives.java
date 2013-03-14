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

import com.google.java.contract.Requires;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import java.util.*;

/**
 * Takes alleles from a variants file and breaks them up (if possible) into more basic/primitive alleles.
 *
 * <p>
 * For now this tool modifies only multi-nucleotide polymorphisms (MNPs) and leaves SNPs, indels, and complex substitutions as is,
 * although one day it may be extended to handle the complex substitution case.
 *
 * This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).
 *
 * Note that this tool modifies only bi-allelic variants.
 *
 * <h2>Input</h2>
 * <p>
 * A variant set with any type of alleles.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A VCF with alleles broken into primitive types.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantsToAllelicPrimitives \
 *   --variant input.vcf \
 *   -o output.vcf
 * </pre>
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
