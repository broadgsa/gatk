package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;

import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 17, 2009
 * Time: 2:48:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProportionOfNonrefBasesSupportingSNP implements VariantAnnotation{
    private String KEY_NAME = "prop_nonref_that_are_snp";

    public String getKeyName() { return KEY_NAME; }

    public VCFInfoHeaderLine getDescription() {
        return new VCFInfoHeaderLine(KEY_NAME,
                        1,VCFInfoHeaderLine.INFO_TYPE.Float,"Simple proportion of non-reference bases that are the SNP base");
    }

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, Variation var) {
        if ( ! var.isSNP() || ! var.isBiallelic() ) {
            return null;
        } else {
            Pair<Integer,Integer> totalNonref_totalSNP = new Pair<Integer,Integer>(0,0);
            for ( String sample : context.keySet() ) {
                ReadBackedPileup pileup = context.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getPileup();
                totalNonref_totalSNP = getNonrefAndSNP(pileup, ref.getBase(), var.getAlternativeBaseForSNP(), totalNonref_totalSNP);

            }
            double p = getProportionOfNonrefBasesThatAreSNP(totalNonref_totalSNP);
            return String.format("%f", p );
        }
    }

    private Pair<Integer,Integer> getNonrefAndSNP(ReadBackedPileup p, char ref, char snp, Pair<Integer,Integer> totals) {
        int[] counts = p.getBaseCounts();
        int nonrefCounts = 0;
        int snpCounts = counts[BaseUtils.simpleBaseToBaseIndex(snp)];
        for ( char c : BaseUtils.BASES ) {
            if ( ! BaseUtils.basesAreEqual((byte) c, (byte) ref) ) {
                nonrefCounts += counts[BaseUtils.simpleBaseToBaseIndex(c)];
            }
        }

        totals.first+=nonrefCounts;
        totals.second+=snpCounts;
        return totals;
    }

    private double getProportionOfNonrefBasesThatAreSNP( Pair<Integer,Integer> totalNonref_totalSNP ) {
        return ( 1.0 + totalNonref_totalSNP.second ) / (1.0 + totalNonref_totalSNP.first );
    }
}
