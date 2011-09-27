package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
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
 * Fraction of all reads across samples that have mapping quality zero
 */
public class MappingQualityZeroFraction extends InfoFieldAnnotation implements ExperimentalAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int mq0 = 0;
        int depth = 0;
        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            AlignmentContext context = sample.getValue();
            depth += context.size();
             ReadBackedPileup pileup = null;
            if (context.hasExtendedEventPileup())
                pileup = context.getExtendedEventPileup();
            else if (context.hasBasePileup())
                pileup = context.getBasePileup();

            if (pileup != null) {
                for (PileupElement p : pileup ) {
                    if ( p.getMappingQual() == 0 )
                        mq0++;
                }
            }
        }
        if (depth > 0) {
            double mq0f = (double)mq0 / (double )depth;

            Map<String, Object> map = new HashMap<String, Object>();
            map.put(getKeyNames().get(0), String.format("%1.4f", mq0f));
            return map;
        }
        else
            return null;
    }

    public List<String> getKeyNames() { return Arrays.asList("MQ0Fraction"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Integer, "Fraction of Mapping Quality Zero Reads")); }
}