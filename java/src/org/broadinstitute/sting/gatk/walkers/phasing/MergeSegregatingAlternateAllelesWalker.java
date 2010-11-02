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

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

import static org.broadinstitute.sting.utils.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci, and merges consecutive sites if they segregate in all samples in the ROD.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE}, referenceMetaData = @RMD(name = "variant", type = ReferenceOrderedDatum.class))
@By(DataSource.REFERENCE_ORDERED_DATA)

public class MergeSegregatingAlternateAllelesWalker extends RodWalker<Integer, Integer> {

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter writer = null;
    private MergePhasedSegregatingAlternateAllelesVCFWriter vcMergerWriter = null;

    @Argument(fullName = "maxGenomicDistanceForMNP", shortName = "maxDistMNP", doc = "The maximum reference-genome distance between consecutive heterozygous sites to permit merging phased VCF records into a MNP record; [default:1]", required = false)
    protected int maxGenomicDistanceForMNP = 1;

    @Argument(fullName = "useSingleSample", shortName = "useSample", doc = "Only output genotypes for the single sample given; [default:use all samples]", required = false)
    protected String useSingleSample = null;

    @Hidden
    @Argument(fullName = "emitOnlyMergedRecords", shortName = "emitOnlyMerged", doc = "Only output records that resulted from merging [For DEBUGGING purposes only - DO NOT USE, since it disregards the semantics of '|' as 'phased relative to previous non-filtered VC']; [default:false]", required = false)
    protected boolean emitOnlyMergedRecords = false;

    @Argument(fullName = "disablePrintAltAlleleStats", shortName = "noAlleleStats", doc = "Should the print-out of alternate allele statistics be disabled?; [default:false]", required = false)
    protected boolean disablePrintAlternateAlleleStatistics = false;

    private LinkedList<String> rodNames = null;

    public void initialize() {
        rodNames = new LinkedList<String>();
        rodNames.add("variant");

        initializeVcfWriter();
    }

    private void initializeVcfWriter() {
        // false <-> don't take control of writer, since didn't create it:
        vcMergerWriter = new MergePhasedSegregatingAlternateAllelesVCFWriter(writer, getToolkit().getArguments().referenceFile, maxGenomicDistanceForMNP, useSingleSample, emitOnlyMergedRecords, logger, false, !disablePrintAlternateAlleleStatistics);
        writer = null; // so it can't be accessed directly [i.e., not through vcMergerWriter]

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), rodNames);
        vcMergerWriter.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(rodNameToHeader.get(rodNames.get(0)).getGenotypeSamples())));
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

        boolean requireStartHere = true; // only see each VariantContext once
        boolean takeFirstOnly = false; // take as many entries as the VCF file has
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly))
            writeVCF(vc);

        return 0;
    }

    private void writeVCF(VariantContext vc) {
        byte refBase;
        if (!vc.isIndel()) {
            Allele varAllele = vc.getReference();
            refBase = SNPallelePair.getSingleBase(varAllele);
        }
        else {
            refBase = vc.getReferenceBaseForIndel();
        }

        vcMergerWriter.add(vc, refBase);
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

        if (useSingleSample != null)
            System.out.println("Only considered single sample: " + useSingleSample);

        System.out.println("Number of successive pairs of records (any distance): " + vcMergerWriter.getNumRecordsAttemptToMerge());
        System.out.println("Number of potentially merged records (distance <= " + maxGenomicDistanceForMNP + "): " + vcMergerWriter.getNumRecordsWithinDistance());
        System.out.println("Number of records merged [all samples are mergeable, some sample has a MNP of ALT alleles]: " + vcMergerWriter.getNumMergedRecords());
        System.out.println(vcMergerWriter.getAltAlleleStats());
    }
}