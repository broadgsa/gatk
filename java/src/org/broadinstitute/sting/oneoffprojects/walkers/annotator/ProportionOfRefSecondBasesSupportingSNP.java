package org.broadinstitute.sting.oneoffprojects.walkers.annotator;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;

import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 17, 2009
 * Time: 2:18:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProportionOfRefSecondBasesSupportingSNP implements VariantAnnotation {
    private String KEY_NAME = "ref_2bb_snp_prop";
    private boolean USE_MAPQ0_READS = false;

    public String getKeyName() { return KEY_NAME; }

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> context, VariantContext vc) {
        if ( ! vc.isSNP() || ! vc.isBiallelic() ) {
            return null;
        } else {
            Pair<Integer,Integer> totalAndSNPSupporting = new Pair<Integer,Integer>(0,0);
            for ( String sample : context.keySet() ) {
                ReadBackedPileup pileup = context.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
                totalAndSNPSupporting = getTotalRefAndSNPSupportCounts(pileup, ref.getBase(), vc.getAlternateAllele(0).toString().charAt(0), totalAndSNPSupporting);

            }
            if ( totalAndSNPSupporting.equals(new Pair<Integer,Integer>(0,0)) )
                return null;
            double p = getProportionOfRefSecondaryBasesSupportingSNP(totalAndSNPSupporting);
            return String.format("%f", p );
        }
    }

    private double getProportionOfRefSecondaryBasesSupportingSNP(Pair<Integer,Integer> tRef_snpSupport) {
        return ( 1.0 + tRef_snpSupport.second) / (1.0 + tRef_snpSupport.first );
    }

    private Pair<Integer,Integer> getTotalRefAndSNPSupportCounts(ReadBackedPileup p, char ref, char snp, Pair<Integer,Integer> refAndSNPCounts) {
        int nRefBases = 0;
        int nSecondBasesSupportingSNP = 0;
        for (PileupElement e : p ) {
            if ( BaseUtils.basesAreEqual( e.getBase(), (byte) ref ) ) {
                if ( BaseUtils.isRegularBase(e.getSecondBase()) ) {
                    nRefBases++;
                    if ( BaseUtils.basesAreEqual( e.getSecondBase(), (byte) snp ) ) {
                        nSecondBasesSupportingSNP++;
                    }
                }
            }
        }

        refAndSNPCounts.first+=nRefBases;
        refAndSNPCounts.second+=nSecondBasesSupportingSNP;
        return refAndSNPCounts;
    }

    public VCFInfoHeaderLine getDescription() {
        return new VCFInfoHeaderLine(KEY_NAME,
                        1,VCFInfoHeaderLine.INFO_TYPE.Float,"Simple proportion of second best base calls for reference base that support the SNP base");
    }
}
