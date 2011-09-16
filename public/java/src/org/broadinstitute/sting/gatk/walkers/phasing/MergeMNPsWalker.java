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

package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

import static org.broadinstitute.sting.utils.codecs.vcf.VCFUtils.getVCFHeadersFromRods;


/**
 * Walks along all variant ROD loci, and merges consecutive sites if they segregate in all samples in the ROD.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})
@By(DataSource.REFERENCE_ORDERED_DATA)

public class MergeMNPsWalker extends RodWalker<Integer, Integer> {

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter writer = null;
    private MergeSegregatingAlternateAllelesVCFWriter vcMergerWriter = null;

    @Argument(fullName = "maxGenomicDistanceForMNP", shortName = "maxDistMNP", doc = "The maximum reference-genome distance between consecutive heterozygous sites to permit merging phased VCF records into a MNP record; [default:1]", required = false)
    protected int maxGenomicDistanceForMNP = 1;

    @Input(fullName="variant", shortName = "V", doc="Select variants from this VCF file", required=true)
    public RodBinding<VariantContext> variants;

    public void initialize() {
        initializeVcfWriter();
    }

    private void initializeVcfWriter() {
        // false <-> don't take control of writer, since didn't create it:
        vcMergerWriter = new MergeSegregatingAlternateAllelesVCFWriter(writer, getToolkit().getGenomeLocParser(), getToolkit().getArguments().referenceFile, maxGenomicDistanceForMNP, logger, false);
        writer = null; // so it can't be accessed directly [i.e., not through vcMergerWriter]

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), Arrays.asList(variants.getName()));
        vcMergerWriter.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(rodNameToHeader.get(variants.getName()).getGenotypeSamples())));
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * For each site, send it to be (possibly) merged with previously observed sites.
     *
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return dummy Integer
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        for (VariantContext vc : tracker.getValues(variants, context.getLocation()))
            writeVCF(vc);

        return 0;
    }

    private void writeVCF(VariantContext vc) {
        WriteVCF.writeVCF(vc, vcMergerWriter, logger);
    }

    public Integer reduce(Integer result, Integer total) {
        if (result == null)
            return total;

        return total + result;
    }

    /**
     * Release any VariantContexts not yet processed.
     *
     * @param result Empty for now...
     */
    public void onTraversalDone(Integer result) {
        vcMergerWriter.close();

        System.out.println("Number of successive pairs of records: " + vcMergerWriter.getNumRecordsAttemptToMerge());
        System.out.println("Number of potentially merged records (" + vcMergerWriter.getVcMergeRule() + "): " + vcMergerWriter.getNumRecordsSatisfyingMergeRule());
        System.out.println("Number of records merged ("+ vcMergerWriter.getAlleleMergeRule() + "): " + vcMergerWriter.getNumMergedRecords());
        System.out.println(vcMergerWriter.getAltAlleleStats());
    }
}