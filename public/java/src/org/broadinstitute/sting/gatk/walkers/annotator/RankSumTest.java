package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.utils.MannWhitneyU;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Abstract root for all RankSum based annotations
 */
public abstract class RankSumTest extends InfoFieldAnnotation implements StandardAnnotation {
    static final double INDEL_LIKELIHOOD_THRESH = 0.1;
    static final boolean DEBUG = false;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;
         
        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;


        final ArrayList<Double> refQuals = new ArrayList<Double>();
        final ArrayList<Double> altQuals = new ArrayList<Double>();

        if (vc.isSNP() && vc.isBiallelic()) {
            // todo - no current support for multiallelic snps
            for ( final Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
                final AlignmentContext context = stratifiedContexts.get(genotype.getKey());
                if ( context == null ) {
                    continue;
                }
                fillQualsFromPileup(ref.getBase(), vc.getAlternateAllele(0).getBases()[0], context.getBasePileup(), refQuals, altQuals);
            }
        }
        else if (vc.isIndel() || vc.isMixed()) {

            for ( final Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
                final AlignmentContext context = stratifiedContexts.get(genotype.getKey());
                if ( context == null ) {
                    continue;
                }

                ReadBackedPileup pileup = null;
                if (context.hasExtendedEventPileup())
                    pileup = context.getExtendedEventPileup();
                else if (context.hasBasePileup())
                    pileup = context.getBasePileup();

                if (pileup == null)
                    continue;

                if (IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap() == null ||
                        IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap().size() == 0)
                    return null;

                fillIndelQualsFromPileup(pileup, refQuals, altQuals);
            }
        }
        else
            return null;

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();
        for ( final Double qual : altQuals ) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET1);
        }
        for ( final Double qual : refQuals ) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET2);
        }

        if (DEBUG) {
            System.out.format("%s, REF QUALS:",this.getClass().getName());
            for ( final Double qual : refQuals )
                System.out.format("%4.1f ",qual);
            System.out.println();
            System.out.format("%s, ALT QUALS:",this.getClass().getName());            
            for ( final Double qual : altQuals )
                System.out.format("%4.1f ",qual);
            System.out.println();

        }
        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final Pair<Double,Double> testResults = mannWhitneyU.runOneSidedTest( MannWhitneyU.USet.SET1 );

        final Map<String, Object> map = new HashMap<String, Object>();
        if ( ! Double.isNaN(testResults.first) )
            map.put(getKeyNames().get(0), String.format("%.3f", testResults.first));
        return map;

    }

    protected abstract void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals);
    protected abstract void fillIndelQualsFromPileup(ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals);

    protected static boolean isUsableBase( final PileupElement p ) {
        return !( p.isDeletion() ||
                  p.getMappingQual() == 0 ||
                  p.getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
                  ((int)p.getQual()) < QualityUtils.MIN_USABLE_Q_SCORE ); // need the unBAQed quality score here
    }
}