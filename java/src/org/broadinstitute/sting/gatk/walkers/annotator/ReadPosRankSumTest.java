package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/30/11
 */

public class ReadPosRankSumTest extends RankSumTest {

    public List<String> getKeyNames() { return Arrays.asList("ReadPosRankSum"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ReadPosRankSum", 1, VCFHeaderLineType.Float, "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias")); }

    protected void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        for ( final PileupElement p : pileup ) {
            if( isUsableBase(p) ) {
                int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p.getOffset(), 0, 0);
                final int numAlignedBases = AlignmentUtils.getNumAlignedBases(p.getRead());
                if( readPos > numAlignedBases / 2 ) {
                    readPos = numAlignedBases - ( readPos + 1 );
                }

                if ( p.getBase() == ref ) {
                    refQuals.add( (double)readPos );
                } else if ( p.getBase() == alt ) {
                    altQuals.add( (double)readPos );
                }
            }
        }
    }
    protected void fillIndelQualsFromPileup(ReadBackedPileup pileup, List<Double> refQuals, List<Double> altQuals) {
        // equivalent is whether indel likelihoods for reads corresponding to ref allele are more likely than reads corresponding to alt allele
        // to classify a pileup element as ref or alt, we look at the likelihood associated with the allele associated to this element.
        // A pileup element then has a list of pairs of form (Allele, likelihood of this allele).
        // To classify a pileup element as Ref or Alt, we look at the likelihood of corresponding alleles.
        // If likelihood of ref allele > highest likelihood of all alt alleles  + epsilon, then this pielup element is "ref"
        // otherwise  if highest alt allele likelihood is > ref likelihood + epsilon, then this pileup element it "alt"
        final HashMap<PileupElement,LinkedHashMap<Allele,Double>> indelLikelihoodMap = IndelGenotypeLikelihoodsCalculationModel.getIndelLikelihoodMap();
        for (final PileupElement p: pileup) {
            if (indelLikelihoodMap.containsKey(p)) {
                // retrieve likelihood information corresponding to this read
                LinkedHashMap<Allele,Double> el = indelLikelihoodMap.get(p);
                // by design, first element in LinkedHashMap was ref allele
                double refLikelihood=0.0, altLikelihood=Double.NEGATIVE_INFINITY;

                for (Allele a : el.keySet()) {

                    if (a.isReference())
                        refLikelihood =el.get(a);
                    else {
                        double like = el.get(a);
                        if (like >= altLikelihood)
                            altLikelihood = like;
                    }
                }
                int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p.getOffset(), 0, 0);
                final int numAlignedBases = AlignmentUtils.getNumAlignedBases(p.getRead());
                if( readPos > numAlignedBases / 2 ) {
                    readPos = numAlignedBases - ( readPos + 1 );
                }

                if (refLikelihood > (altLikelihood + INDEL_LIKELIHOOD_THRESH))
                    refQuals.add((double)readPos);
                else if (altLikelihood > (refLikelihood + INDEL_LIKELIHOOD_THRESH))
                    altQuals.add((double)readPos);


            }
        }
    }

}
