package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
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
public abstract class RankSumTest extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {
    static final boolean DEBUG = false;

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        // either stratifiedContexts or  stratifiedPerReadAlleleLikelihoodMap has to be non-null

        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || genotypes.size() == 0)
            return null;

        final ArrayList<Double> refQuals = new ArrayList<Double>();
        final ArrayList<Double> altQuals = new ArrayList<Double>();

        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            PerReadAlleleLikelihoodMap indelLikelihoodMap = null;
            ReadBackedPileup pileup = null;
            if (stratifiedPerReadAlleleLikelihoodMap != null && !stratifiedPerReadAlleleLikelihoodMap.isEmpty()) {
                indelLikelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(genotype.getSampleName());
                if (indelLikelihoodMap == null)
                    continue;
                if (indelLikelihoodMap.isEmpty())
                    continue;
            }
            else if (stratifiedContexts != null) {
                final AlignmentContext context = stratifiedContexts.get(genotype.getSampleName());
                if ( context == null )
                    continue;
                pileup = context.getBasePileup();
            }
            fillQualsFromPileup(vc.getAlleles(), vc.getStart(), pileup, indelLikelihoodMap, refQuals, altQuals );
        }
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

    protected abstract void fillQualsFromPileup(final List<Allele> alleles,
                                                final int refLoc,
                                                final ReadBackedPileup readBackedPileup,
                                                final PerReadAlleleLikelihoodMap alleleLikelihoodMap,
                                                final List<Double> refQuals,
                                                final List<Double> altQuals);

    /**
     * Can the base in this pileup element be used in comparative tests between ref / alt bases?
     *
     * Note that this function by default does not allow deletion pileup elements
     *
     * @param p the pileup element to consider
     * @return true if this base is part of a meaningful read for comparison, false otherwise
     */
    public static boolean isUsableBase(final PileupElement p) {
        return isUsableBase(p, false);
    }

    /**
     * Can the base in this pileup element be used in comparative tests between ref / alt bases?
     *
     * @param p the pileup element to consider
     * @param allowDeletions if true, allow p to be a deletion base
     * @return true if this base is part of a meaningful read for comparison, false otherwise
     */
    public static boolean isUsableBase(final PileupElement p, final boolean allowDeletions) {
        return !(p.isInsertionAtBeginningOfRead() ||
                 (! allowDeletions && p.isDeletion()) ||
                 p.getMappingQual() == 0 ||
                 p.getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
                 ((int) p.getQual()) < QualityUtils.MIN_USABLE_Q_SCORE); // need the unBAQed quality score here
    }
}