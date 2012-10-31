package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;

import java.util.*;


/**
 * The u-based z-approximation from the Mann-Whitney Rank Sum Test for base qualities (ref bases vs. bases of the alternate allele).
 * Note that the base quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.
 */
public class BaseQualityRankSumTest extends RankSumTest implements StandardAnnotation {
    public List<String> getKeyNames() { return Arrays.asList("BaseQRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("BaseQRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities")); }

    protected void fillQualsFromPileup(final List<Allele> allAlleles, final int refLoc,
                                       final ReadBackedPileup pileup,
                                       final PerReadAlleleLikelihoodMap alleleLikelihoodMap,
                                       final List<Double> refQuals, final List<Double> altQuals){

        if (alleleLikelihoodMap == null) {
            // use fast SNP-based version if we don't have per-read allele likelihoods
            for ( final PileupElement p : pileup ) {
                if ( isUsableBase(p) ) {
                    if ( allAlleles.get(0).equals(Allele.create(p.getBase(),true)) ) {
                        refQuals.add((double)p.getQual());
                    } else if ( allAlleles.contains(Allele.create(p.getBase()))) {
                        altQuals.add((double)p.getQual());
                    }
                }
            }
            return;
        }

        for (Map<Allele,Double> el : alleleLikelihoodMap.getLikelihoodMapValues()) {
            final Allele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el);
            if (a.isNoCall())
                continue; // read is non-informative
            if (a.isReference())
                refQuals.add(-10.0*(double)el.get(a));
            else if (allAlleles.contains(a))
                altQuals.add(-10.0*(double)el.get(a));


        }
    }


}