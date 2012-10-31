package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFStandardHeaderLines;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Total count across all samples of mapping quality zero reads
 */
public class MappingQualityZero extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ((vc.isSNP() || !vc.isVariant()) && stratifiedContexts != null)
            return annotatePileup(ref, stratifiedContexts, vc);
        else if (stratifiedPerReadAlleleLikelihoodMap != null && vc.isVariant())
            return annotateWithLikelihoods(stratifiedPerReadAlleleLikelihoodMap, vc);
        else
            return null;
    }

    private Map<String, Object> annotatePileup(final ReferenceContext ref,
                                               final Map<String, AlignmentContext> stratifiedContexts,
                                               final VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int mq0 = 0;
        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            final AlignmentContext context = sample.getValue();
            final ReadBackedPileup pileup = context.getBasePileup();
            for (PileupElement p : pileup ) {
                if ( p.getMappingQual() == 0 )
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", mq0));
        return map;
    }

    private Map<String, Object> annotateWithLikelihoods(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                        final VariantContext vc) {
        if (stratifiedPerReadAlleleLikelihoodMap == null)
            return null;

        int mq0 = 0;
        for ( PerReadAlleleLikelihoodMap likelihoodMap : stratifiedPerReadAlleleLikelihoodMap.values() ) {
            for (GATKSAMRecord read : likelihoodMap.getLikelihoodReadMap().keySet()) {

                 if (read.getMappingQuality() == 0 )
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", mq0));
        return map;
    }


    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.MAPPING_QUALITY_ZERO_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}