package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Triplet annotation: fraction of MAQP == 0, MAPQ < 10, and count of all mapped reads
 */
public class LowMQ extends InfoFieldAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        double mq0 = 0;
		double mq10 = 0;
		double total = 0;
        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() )
		{
            ReadBackedPileup pileup = sample.getValue().getBasePileup();
            for (PileupElement p : pileup ) 
			{
                if ( p.getMappingQual() == 0 )  { mq0 += 1; }
                if ( p.getMappingQual() <= 10 ) { mq10 += 1; }
				total += 1; 
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.04f,%.04f,%.00f", mq0/total, mq10/total, total));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("LowMQ"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 3, VCFHeaderLineType.Float, "3-tuple: <fraction of reads with MQ=0>,<fraction of reads with MQ<=10>,<total number of reads>")); }
}
