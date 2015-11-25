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

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import htsjdk.samtools.util.IntervalTree;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

@Analysis(description = "1000 Genomes Phase I summary of variants table")
public class VariantSummary extends VariantEvaluator implements StandardEval {
    final protected static Logger logger = Logger.getLogger(VariantSummary.class);

    /** Indels with size greater than this value are tallied in the CNV column */
    private final static int MAX_INDEL_LENGTH = 50;
    private final static double MIN_CNV_OVERLAP = 0.5;

    public enum Type {
        SNP, INDEL, CNV
    }

    Map<String, IntervalTree<GenomeLoc>> knownCNVs = null;

    // basic counts on various rates found
    @DataPoint(description = "Number of samples", format = "%d")
    public long nSamples = 0;

    @DataPoint(description = "Number of processed loci", format = "%d")
    public long nProcessedLoci = 0;

    @DataPoint(description = "Number of SNPs", format = "%d")
    public long nSNPs = 0;
    @DataPoint(description = "Overall TiTv ratio", format = "%.2f")
    public double TiTvRatio = 0;
    @DataPoint(description = "SNP Novelty Rate", format = "%s")
    public String SNPNoveltyRate = "NA";
    @DataPoint(description = "Mean number of SNPs per individual", format = "%d")
    public long nSNPsPerSample = 0;
    @DataPoint(description = "Mean TiTv ratio per individual", format = "%.2f")
    public double TiTvRatioPerSample = 0;
    @DataPoint(description = "Mean depth of coverage per sample at SNPs", format = "%.1f")
    public double SNPDPPerSample = 0;

    @DataPoint(description = "Number of Indels", format = "%d")
    public long nIndels = 0;
    @DataPoint(description = "Indel Novelty Rate", format = "%s")
    public String IndelNoveltyRate = "NA";
    @DataPoint(description = "Mean number of Indels per individual", format = "%d")
    public long nIndelsPerSample = 0;
    @DataPoint(description = "Mean depth of coverage per sample at Indels", format = "%.1f")
    public double IndelDPPerSample = 0;

    @DataPoint(description = "Number of SVs", format = "%d")
    public long nSVs = 0;
    @DataPoint(description = "SV Novelty Rate", format = "%s")
    public String SVNoveltyRate = "NA";
    @DataPoint(description = "Mean number of SVs per individual", format = "%d")
    public long nSVsPerSample = 0;

    TypeSampleMap allVariantCounts, knownVariantCounts;
    TypeSampleMap countsPerSample;
    TypeSampleMap transitionsPerSample, transversionsPerSample;
    TypeSampleMap depthPerSample;

    private final static String ALL = "ALL";

    private class TypeSampleMap extends EnumMap<Type, Map<String, Integer>> {
        public TypeSampleMap(final Collection<String> samples) {
            super(Type.class);
            for ( Type type : Type.values() ) {
                Map<String, Integer> bySample = new HashMap<String, Integer>(samples.size());
                for ( final String sample : samples ) {
                    bySample.put(sample, 0);
                }
                bySample.put(ALL, 0);
                this.put(type, bySample);
            }
        }

        public final void inc(final Type type, final String sample) {
            final int count = this.get(type).get(sample);
            get(type).put(sample, count + 1);
        }

        public final int all(Type type) {
            return get(type).get(ALL);
        }

        public final int meanValue(Type type) {
            long sum = 0;
            int n = 0;
            for ( final Map.Entry<String, Integer> pair : get(type).entrySet() ) {
                if ( pair.getKey() != ALL)  { // truly must be string ==
                    n++;
                    sum += pair.getValue();
                }
            }
            return (int)(Math.round(sum / (1.0 * n)));
        }

        public final double ratioValue(Type type, TypeSampleMap denoms, boolean allP) {
            double sum = 0;
            int n = 0;
            for ( final String sample : get(type).keySet() ) {
                if ( (allP && sample == ALL) || (!allP && sample != ALL) ) { // truly must be string ==
                    final long num = get(type).get(sample);
                    final long denom = denoms.get(type).get(sample);
                    sum += ratio(num, denom);
                    n++;
                }
            }

            return n > 0 ? sum / (1.0 * n) : 0.0;
        }
    }


    public void initialize(VariantEval walker) {
        super.initialize(walker);

        nSamples = walker.getSampleNamesForEvaluation().size();
        countsPerSample = new TypeSampleMap(walker.getSampleNamesForEvaluation());
        transitionsPerSample = new TypeSampleMap(walker.getSampleNamesForEvaluation());
        transversionsPerSample = new TypeSampleMap(walker.getSampleNamesForEvaluation());
        allVariantCounts = new TypeSampleMap(walker.getSampleNamesForEvaluation());
        knownVariantCounts = new TypeSampleMap(walker.getSampleNamesForEvaluation());
        depthPerSample = new TypeSampleMap(walker.getSampleNamesForEvaluation());

        if ( walker.knownCNVsFile != null ) {
            knownCNVs = walker.createIntervalTreeByContig(walker.knownCNVsFile);
            final List<GenomeLoc> locs = walker.knownCNVsFile.getIntervals(walker.getToolkit().getGenomeLocParser());
            logger.info(String.format("Creating known CNV list %s containing %d intervals covering %d bp",
                    walker.knownCNVsFile.getSource(), locs.size(), IntervalUtils.intervalSize(locs)));
        }
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    private Type getType(VariantContext vc) {
        switch (vc.getType()) {
            case SNP:
                return Type.SNP;
            case INDEL:
                for ( int l : vc.getIndelLengths() )
                    if ( Math.abs(l) > MAX_INDEL_LENGTH )
                        return Type.CNV;
                return Type.INDEL;
            case SYMBOLIC:
                return Type.CNV;
            default:
                //throw new UserException.BadInput("Unexpected variant context type: " + vc);
                return null;
        }
    }

    private boolean overlapsKnownCNV(VariantContext cnv) {
        if ( knownCNVs != null ) {
            final GenomeLoc loc = getWalker().getToolkit().getGenomeLocParser().createGenomeLoc(cnv);
            IntervalTree<GenomeLoc> intervalTree = knownCNVs.get(loc.getContig());

            final Iterator<IntervalTree.Node<GenomeLoc>> nodeIt = intervalTree.overlappers(loc.getStart(), loc.getStop());
            while ( nodeIt.hasNext() ) {
                final double overlapP = loc.reciprocialOverlapFraction(nodeIt.next().getValue());
                if ( overlapP > MIN_CNV_OVERLAP )
                    return true;
            }
        }

        return false;
    }

    public void update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval == null || (getWalker().ignoreAC0Sites() && eval.isMonomorphicInSamples()) )
            return;

        final Type type = getType(eval);
        if ( type == null )
            return;

        TypeSampleMap titvTable = null;

        // update DP, if possible
        if ( eval.hasAttribute(VCFConstants.DEPTH_KEY) )
            depthPerSample.inc(type, ALL);

        // update counts
        allVariantCounts.inc(type, ALL);

        // type specific calculations
        if ( type == Type.SNP && eval.isBiallelic() ) {
            titvTable = GATKVariantContextUtils.isTransition(eval) ? transitionsPerSample : transversionsPerSample;
            titvTable.inc(type, ALL);
        }

        // novelty calculation
        if ( comp != null || (type == Type.CNV && overlapsKnownCNV(eval)))
            knownVariantCounts.inc(type, ALL);

        // per sample metrics
        for (final Genotype g : eval.getGenotypes()) {
            if ( ! g.isNoCall() && ! g.isHomRef() ) {
                countsPerSample.inc(type, g.getSampleName());

                // update transition / transversion ratio
                if ( titvTable != null ) titvTable.inc(type, g.getSampleName());

                if ( g.hasDP() )
                    depthPerSample.inc(type, g.getSampleName());
            }
        }
    }

    private String noveltyRate(Type type) {
        final int all = allVariantCounts.all(type);
        final int known = knownVariantCounts.all(type);
        return Utils.formattedNoveltyRate(known, all);
    }

    public void finalizeEvaluation() {
        nProcessedLoci = getWalker().getnProcessedLoci();
        nSNPs = allVariantCounts.all(Type.SNP);
        nIndels = allVariantCounts.all(Type.INDEL);
        nSVs = allVariantCounts.all(Type.CNV);

        TiTvRatio = transitionsPerSample.ratioValue(Type.SNP, transversionsPerSample, true);
        TiTvRatioPerSample = transitionsPerSample.ratioValue(Type.SNP, transversionsPerSample, false);

        nSNPsPerSample = countsPerSample.meanValue(Type.SNP);
        nIndelsPerSample = countsPerSample.meanValue(Type.INDEL);
        nSVsPerSample = countsPerSample.meanValue(Type.CNV);

        SNPNoveltyRate = noveltyRate(Type.SNP);
        IndelNoveltyRate = noveltyRate(Type.INDEL);
        SVNoveltyRate = noveltyRate(Type.CNV);

        SNPDPPerSample = depthPerSample.meanValue(Type.SNP);
        IndelDPPerSample = depthPerSample.meanValue(Type.INDEL);
    }
}