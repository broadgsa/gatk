package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;


public class LowMQ implements VariantAnnotation {

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        double mq0 = 0;
		double mq10 = 0;
		double total = 0;
        for ( String sample : stratifiedContexts.keySet() ) 
		{
            ReadBackedPileup pileup = stratifiedContexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for (PileupElement p : pileup ) 
			{
                if ( p.getMappingQual() == 0 )  { mq0 += 1; }
                if ( p.getMappingQual() <= 10 ) { mq10 += 1; }
				total += 1; 
            }
        }
        return String.format("%.04f,%.04f,%.00f", mq0/total, mq10/total, total);
    }

    public String getKeyName() { return "LowMQ"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(getKeyName(), 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "3-tuple: <fraction of reads with MQ=0>,<fraction of reads with MQ<=10>,<total nubmer of reads>"); }
}
