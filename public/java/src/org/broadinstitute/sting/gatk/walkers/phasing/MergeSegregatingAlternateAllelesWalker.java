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

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

import static org.broadinstitute.sting.utils.codecs.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci, and merges consecutive sites if some sample has segregating alt alleles in the ROD.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})
@By(DataSource.REFERENCE_ORDERED_DATA)

public class MergeSegregatingAlternateAllelesWalker extends RodWalker<Integer, Integer> {

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter writer = null;
    private MergeSegregatingAlternateAllelesVCFWriter vcMergerWriter = null;

    @Argument(fullName = "maxGenomicDistance", shortName = "maxDist", doc = "The maximum reference-genome distance between consecutive heterozygous sites to permit merging phased VCF records; [default:1]", required = false)
    protected int maxGenomicDistance = 1;

    @Argument(fullName = "useSingleSample", shortName = "useSample", doc = "Only output genotypes for the single sample given; [default:use all samples]", required = false)
    protected String useSingleSample = null;

    @Hidden
    @Argument(fullName = "emitOnlyMergedRecords", shortName = "emitOnlyMerged", doc = "Only output records that resulted from merging [For DEBUGGING purposes only - DO NOT USE, since it disregards the semantics of '|' as 'phased relative to previous non-filtered VC']; [default:false]", required = false)
    protected boolean emitOnlyMergedRecords = false;

    @Argument(fullName = "disablePrintAltAlleleStats", shortName = "noAlleleStats", doc = "Should the print-out of alternate allele statistics be disabled?; [default:false]", required = false)
    protected boolean disablePrintAlternateAlleleStatistics = false;

    public final static String IGNORE_REFSEQ = "IGNORE";
    public final static String UNION_REFSEQ = "UNION";
    public final static String INTERSECT_REFSEQ = "INTERSECT";

    @Argument(fullName = "mergeBasedOnRefSeqAnnotation", shortName = "mergeBasedOnRefSeqAnnotation", doc = "'Should merging be performed if two sites lie on the same RefSeq sequence in the INFO field {" + IGNORE_REFSEQ + ", " + UNION_REFSEQ + ", " + INTERSECT_REFSEQ + "}; [default:" + IGNORE_REFSEQ + "]", required = false)
    protected String mergeBasedOnRefSeqAnnotation = IGNORE_REFSEQ;

    @Argument(fullName = "dontRequireSomeSampleHasDoubleAltAllele", shortName = "dontRequireSomeSampleHasDoubleAltAllele", doc = "Should the requirement, that SUCCESSIVE records to be merged have at least one sample with a double alternate allele, be relaxed?; [default:false]", required = false)
    protected boolean dontRequireSomeSampleHasDoubleAltAllele = false;

    @Input(fullName="variant", shortName = "V", doc="Select variants from this VCF file", required=true)
    public RodBinding<VariantContext> variants;

    public void initialize() {
        initializeVcfWriter();
    }

    private void initializeVcfWriter() {
        GenomeLocParser genomeLocParser = getToolkit().getGenomeLocParser();

        VariantContextMergeRule vcMergeRule;
        if (mergeBasedOnRefSeqAnnotation.equals(IGNORE_REFSEQ))
            vcMergeRule = new DistanceMergeRule(maxGenomicDistance, genomeLocParser);
        else
            vcMergeRule = new SameGenePlusWithinDistanceMergeRule(maxGenomicDistance, genomeLocParser, mergeBasedOnRefSeqAnnotation);

        VariantContextUtils.AlleleMergeRule alleleMergeRule;
        if (dontRequireSomeSampleHasDoubleAltAllele) // if a pair of VariantContext passes the vcMergeRule, then always merge them if there is a trailing prefix of polymorphisms (i.e., upstream polymorphic site):
            alleleMergeRule = new PrefixPolymorphismMergeAllelesRule();
        else
            alleleMergeRule = new ExistsDoubleAltAlleleMergeRule();

        // false <-> don't take control of writer, since didn't create it:
        vcMergerWriter = new MergeSegregatingAlternateAllelesVCFWriter(writer, genomeLocParser, getToolkit().getArguments().referenceFile, vcMergeRule, alleleMergeRule, useSingleSample, emitOnlyMergedRecords, logger, false, !disablePrintAlternateAlleleStatistics);
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

        if (useSingleSample != null)
            System.out.println("Only considered single sample: " + useSingleSample);

        System.out.println("Number of successive pairs of records: " + vcMergerWriter.getNumRecordsAttemptToMerge());
        System.out.println("Number of potentially merged records (" + vcMergerWriter.getVcMergeRule() + "): " + vcMergerWriter.getNumRecordsSatisfyingMergeRule());
        System.out.println("Number of records merged ("+ vcMergerWriter.getAlleleMergeRule() + "): " + vcMergerWriter.getNumMergedRecords());
        System.out.println(vcMergerWriter.getAltAlleleStats());
    }
}


enum MergeBasedOnRefSeqAnnotation {
    UNION_WITH_DIST, INTERSECT_WITH_DIST
}

class SameGenePlusWithinDistanceMergeRule extends DistanceMergeRule {
    private MergeBasedOnRefSeqAnnotation mergeBasedOnRefSeqAnnotation;

    public SameGenePlusWithinDistanceMergeRule(int maxGenomicDistanceForMNP, GenomeLocParser genomeLocParser, String mergeBasedOnRefSeqAnnotation) {
        super(maxGenomicDistanceForMNP, genomeLocParser);

        if (mergeBasedOnRefSeqAnnotation.equals(MergeSegregatingAlternateAllelesWalker.UNION_REFSEQ))
            this.mergeBasedOnRefSeqAnnotation = MergeBasedOnRefSeqAnnotation.UNION_WITH_DIST;
        else if (mergeBasedOnRefSeqAnnotation.equals(MergeSegregatingAlternateAllelesWalker.INTERSECT_REFSEQ))
            this.mergeBasedOnRefSeqAnnotation = MergeBasedOnRefSeqAnnotation.INTERSECT_WITH_DIST;
        else
            throw new UserException("Must provide " + MergeSegregatingAlternateAllelesWalker.IGNORE_REFSEQ + ", " + MergeSegregatingAlternateAllelesWalker.UNION_REFSEQ + ", or " + MergeSegregatingAlternateAllelesWalker.INTERSECT_REFSEQ + " as argument to mergeBasedOnRefSeqAnnotation!");
    }

    public boolean shouldAttemptToMerge(VariantContext vc1, VariantContext vc2) {
        boolean withinDistance = super.shouldAttemptToMerge(vc1, vc2);

        if (mergeBasedOnRefSeqAnnotation == MergeBasedOnRefSeqAnnotation.UNION_WITH_DIST)
            return withinDistance || sameGene(vc1, vc2);
        else // mergeBasedOnRefSeqAnnotation == MergeBasedOnRefSeqAnnotation.INTERSECT_WITH_DIST
            return withinDistance && sameGene(vc1, vc2);
    }

    private boolean sameGene(VariantContext vc1, VariantContext vc2) {
        Set<String> names_vc1 = RefSeqDataParser.getRefSeqNames(vc1);
        Set<String> names_vc2 = RefSeqDataParser.getRefSeqNames(vc2);
        names_vc1.retainAll(names_vc2);

        if (!names_vc1.isEmpty())
            return true;

        // Check refseq.name2:
        Set<String> names2_vc1 = RefSeqDataParser.getRefSeqNames(vc1, true);
        Set<String> names2_vc2 = RefSeqDataParser.getRefSeqNames(vc2, true);
        names2_vc1.retainAll(names2_vc2);

        return !names2_vc1.isEmpty();
    }

    public String toString() {
        return super.toString() + " " + (mergeBasedOnRefSeqAnnotation == MergeBasedOnRefSeqAnnotation.UNION_WITH_DIST ? "OR" : "AND") + " on the same gene";
    }

    public Map<String, Object> addToMergedAttributes(VariantContext vc1, VariantContext vc2) {
        Map<String, Object> addedAttribs = super.addToMergedAttributes(vc1, vc2);
        addedAttribs.putAll(RefSeqDataParser.getMergedRefSeqNameAttributes(vc1, vc2));
        return addedAttribs;
    }
}



class PrefixPolymorphismMergeAllelesRule extends VariantContextUtils.AlleleMergeRule {
    public boolean allelesShouldBeMerged(VariantContext vc1, VariantContext vc2) {
        return vc1.isPolymorphic();
    }

    public String toString() {
        return super.toString() + ", there exists a polymorphism at the start of the merged allele";
    }
}