package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 17, 2009
 * Time: 2:42:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProportionOfSNPSecondBasesSupportingRef implements VariantAnnotation{
    public String KEY_NAME = "SNP_2B_SUPPORT_REF";
    public boolean USE_MAPQ0_READS = false;
    public String debug_file = "/humgen/gsa-scr1/chartl/temporary/ProportionOfRefSecondBasesSupportingSNP.debug.txt";

    public String getKeyName() { return KEY_NAME; }

    public boolean useZeroQualityReads() { return USE_MAPQ0_READS; }

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, Variation var) {
        if ( ! var.isSNP() || ! var.isBiallelic() ) {
            return null;
        } else {
            Pair<Integer,Integer> totalAndSNPSupporting = new Pair<Integer,Integer>(0,0);
            for ( String sample : context.keySet() ) {
                ReadBackedPileup pileup = context.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getPileup();
                totalAndSNPSupporting = getTotalSNPandRefSupporting(pileup, ref.getBase(), var.getAlternativeBaseForSNP(), totalAndSNPSupporting);

            }
            double p = getProportionOfSNPSecondaryBasesSupportingRef(totalAndSNPSupporting);
            return String.format("%f", p );
        }
    }

    public double getProportionOfSNPSecondaryBasesSupportingRef(Pair<Integer,Integer> tSNP_refSupport) {
        return ( 1.0 + tSNP_refSupport.second) / (1.0 + tSNP_refSupport.first );
    }

    public Pair<Integer,Integer> getTotalSNPandRefSupporting(ReadBackedPileup p, char ref, char snp, Pair<Integer,Integer> SNPandRefCounts) {
        int nSNPBases = 0;
        int nSNPBasesSupportingRef = 0;
        for (PileupElement e : p ) {
            if ( BaseUtils.basesAreEqual( e.getBase(), (byte) snp ) ) {
                if ( hasSecondBase(e) ) {
                    nSNPBases++;
                    if ( BaseUtils.basesAreEqual( e.getSecondBase(), (byte) ref ) ) {
                        nSNPBasesSupportingRef++;
                    }
                }
            }
        }

        SNPandRefCounts.first+=nSNPBases;
        SNPandRefCounts.second+=nSNPBasesSupportingRef;
        return SNPandRefCounts;
    }


    public boolean hasSecondBase(PileupElement e) {
        return BaseUtils.isRegularBase(e.getSecondBase());
    }

    public VCFInfoHeaderLine getDescription() {
        return new VCFInfoHeaderLine(KEY_NAME,
                        1,VCFInfoHeaderLine.INFO_TYPE.Float,"Simple proportion of second best base calls for SNP base that support the Ref base");
    }


}
