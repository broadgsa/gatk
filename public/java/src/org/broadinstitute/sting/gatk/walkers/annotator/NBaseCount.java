package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The number of N bases, counting only SOLiD data
 */
public class NBaseCount extends InfoFieldAnnotation {
    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if( stratifiedContexts.size() == 0 )
            return null;

        int countNBaseSolid = 0;
        int countRegularBaseSolid = 0;

        for( final Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            for( final PileupElement p : sample.getValue().getBasePileup()) {
                if( p.getRead().getReadGroup().getPlatform().toUpperCase().contains("SOLID") ) {
                    if( BaseUtils.isNBase( p.getBase() ) ) {
                        countNBaseSolid++;
                    } else if( BaseUtils.isRegularBase( p.getBase() ) ) {
                        countRegularBaseSolid++;
                    }
                }
            }
        }
        final Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.4f", (double)countNBaseSolid / (double)(countNBaseSolid + countRegularBaseSolid + 1)));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("PercentNBaseSolid"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("PercentNBaseSolid", 1, VCFHeaderLineType.Float, "Percentage of N bases in the pileup (counting only SOLiD reads)")); }
}
