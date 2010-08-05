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

package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;

import java.util.*;

/**
 * Extracts subsets of a VCF file like one or more samples, all or only variant loci, all or filtered loci.
 */
public class VariantSubset extends RodWalker<Integer, Integer> {
    @Argument(fullName="sample", shortName="SN", doc="Sample to include (or nothing to specify all samples)", required=false)
    private ArrayList<String> SAMPLES = null;

    @Argument(fullName="includeNonVariants", shortName="INV", doc="Include non-variant loci", required=false)
    private boolean INCLUDE_NON_VARIANTS = false;

    @Argument(fullName="includeFiltered", shortName="IF", doc="Include filtered loci", required=false)
    private boolean INCLUDE_FILTERED = false;

    private VCFWriter writer;

    public void initialize() {
        Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
        metaData.add(new VCFHeaderLine("source", "VariantsToVCF"));
        metaData.add(new VCFHeaderLine("reference", this.getToolkit().getArguments().referenceFile.getAbsolutePath()));

        writer = new VCFWriter(out);
        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList("variant"));

        final VCFHeader vcfHeader = new VCFHeader(metaData, samples);
        writer.writeHeader(vcfHeader);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Collection<VariantContext> VCs = tracker.getAllVariantContexts(ref, null, context.getLocation(), true, false);
        for (VariantContext vc : VCs) {
            VariantContext subset = subsetRecord(vc);

            if ( (vc.isPolymorphic() || INCLUDE_NON_VARIANTS) &&
                    (!subset.isFiltered() || INCLUDE_FILTERED) )
                writer.add(subset, ref.getBase());
        }

        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    private VariantContext subsetRecord(VariantContext vc) {
        if ( SAMPLES == null || SAMPLES.isEmpty() )
            return vc;

        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
        for ( Map.Entry<String, Genotype> genotypePair : vc.getGenotypes().entrySet() ) {
            if ( SAMPLES.contains(genotypePair.getKey()) )
                genotypes.add(genotypePair.getValue());
        }

        return vc.subContextFromGenotypes(genotypes);
    }

    public Integer reduce(Integer sum, Integer value) {
        return 1;
    }

    public void onTraversalDone(Integer sum) {
        writer.close();
    }
}
