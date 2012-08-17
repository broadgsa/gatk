package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 6/28/12
 */

public class ClippingRankSumTest extends RankSumTest {

    public List<String> getKeyNames() { return Arrays.asList("ClippingRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ClippingRankSum", 1, VCFHeaderLineType.Float, "Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases")); }


    protected void fillQualsFromPileup(final List<Allele> allAlleles,
                                       final int refLoc,
                                       final ReadBackedPileup pileup,
                                       final PerReadAlleleLikelihoodMap likelihoodMap, final List<Double> refQuals, final List<Double> altQuals) {
        // todo - only support non-pileup case for now, e.g. active-region based version
        if (pileup != null)
            return;

        for (Map.Entry<PileupElement,Map<Allele,Double>> el : likelihoodMap.getLikelihoodReadMap().entrySet()) {

            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if (a.isNoCall())
                continue; // read is non-informative
            if (a.isReference())
                refQuals.add((double)AlignmentUtils.getNumHardClippedBases(el.getKey().getRead()));
            else if (allAlleles.contains(a))
                altQuals.add((double)AlignmentUtils.getNumHardClippedBases(el.getKey().getRead()));

        }
    }

 }
