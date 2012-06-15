package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.utils.MannWhitneyU;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Abstract root for all RankSum based annotations
 */
public abstract class RankSumTest extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {
    static final double INDEL_LIKELIHOOD_THRESH = 0.1;
    static final boolean DEBUG = false;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if (stratifiedContexts.size() == 0)
            return null;

        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || genotypes.size() == 0)
            return null;

        final ArrayList<Double> refQuals = new ArrayList<Double>();
        final ArrayList<Double> altQuals = new ArrayList<Double>();

        if ( vc.isSNP() ) {
            final List<Byte> altAlleles = new ArrayList<Byte>();
            for ( final Allele a : vc.getAlternateAlleles() )
                altAlleles.add(a.getBases()[0]);

            for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
                final AlignmentContext context = stratifiedContexts.get(genotype.getSampleName());
                if ( context == null )
                    continue;

                fillQualsFromPileup(ref.getBase(), altAlleles, context.getBasePileup(), refQuals, altQuals);
            }
        } else if ( vc.isIndel() || vc.isMixed() ) {

            for (final Genotype genotype : genotypes.iterateInSampleNameOrder()) {
                final AlignmentContext context = stratifiedContexts.get(genotype.getSampleName());
                if (context == null) {
                    continue;
                }

                final ReadBackedPileup pileup = context.getBasePileup();
                if (pileup == null)
                    continue;

                if (IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap() == null ||
                        IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap().size() == 0)
                    return null;

                fillIndelQualsFromPileup(pileup, refQuals, altQuals);
            }
        } else
            return null;

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();
        for (final Double qual : altQuals) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET1);
        }
        for (final Double qual : refQuals) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET2);
        }

        if (DEBUG) {
            System.out.format("%s, REF QUALS:", this.getClass().getName());
            for (final Double qual : refQuals)
                System.out.format("%4.1f ", qual);
            System.out.println();
            System.out.format("%s, ALT QUALS:", this.getClass().getName());
            for (final Double qual : altQuals)
                System.out.format("%4.1f ", qual);
            System.out.println();

        }
        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final Pair<Double, Double> testResults = mannWhitneyU.runOneSidedTest(MannWhitneyU.USet.SET1);

        final Map<String, Object> map = new HashMap<String, Object>();
        if (!Double.isNaN(testResults.first))
            map.put(getKeyNames().get(0), String.format("%.3f", testResults.first));
        return map;
    }

    public Map<String, Object> annotate(Map<String, Map<Allele, List<GATKSAMRecord>>> stratifiedContexts, VariantContext vc) {
        if (stratifiedContexts.size() == 0)
            return null;

        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || genotypes.size() == 0)
            return null;

        final ArrayList<Double> refQuals = new ArrayList<Double>();
        final ArrayList<Double> altQuals = new ArrayList<Double>();

        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            final Map<Allele, List<GATKSAMRecord>> context = stratifiedContexts.get(genotype.getSampleName());
            if ( context == null )
                continue;

            fillQualsFromPileup(vc.getReference(), vc.getAlternateAlleles(), vc.getStart(), context, refQuals, altQuals);
        }

        if ( refQuals.size() == 0 || altQuals.size() == 0 )
            return null;

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();
        for (final Double qual : altQuals) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET1);
        }
        for (final Double qual : refQuals) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET2);
        }

        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final Pair<Double, Double> testResults = mannWhitneyU.runOneSidedTest(MannWhitneyU.USet.SET1);

        final Map<String, Object> map = new HashMap<String, Object>();
        if (!Double.isNaN(testResults.first))
            map.put(getKeyNames().get(0), String.format("%.3f", testResults.first));
        return map;
    }

    protected abstract void fillQualsFromPileup(final Allele ref, final List<Allele> alts, final int refLoc, final Map<Allele, List<GATKSAMRecord>> stratifiedContext, final List<Double> refQuals, List<Double> altQuals);

    protected abstract void fillQualsFromPileup(final byte ref, final List<Byte> alts, final ReadBackedPileup pileup, final List<Double> refQuals, final List<Double> altQuals);

    protected abstract void fillIndelQualsFromPileup(final ReadBackedPileup pileup, final List<Double> refQuals, final List<Double> altQuals);

    protected static boolean isUsableBase(final PileupElement p) {
        return !(p.isInsertionAtBeginningOfRead() ||
                 p.isDeletion() ||
                 p.getMappingQual() == 0 ||
                 p.getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
                 ((int) p.getQual()) < QualityUtils.MIN_USABLE_Q_SCORE); // need the unBAQed quality score here
    }
}