package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFStandardHeaderLines;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;


/**
 * Root Mean Square of the mapping quality of the reads across all samples.
 */
public class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        int totalSize = 0, index = 0;
        int qualities[];
        if (stratifiedContexts != null) {
            if ( stratifiedContexts.size() == 0 )
                return null;

            for ( AlignmentContext context : stratifiedContexts.values() )
                totalSize += context.size();

            qualities = new int[totalSize];

            for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
                AlignmentContext context = sample.getValue();
                for (PileupElement p : context.getBasePileup() )
                    index = fillMappingQualitiesFromPileupAndUpdateIndex(p, index, qualities);
            }
        }
        else if (perReadAlleleLikelihoodMap != null) {
            if ( perReadAlleleLikelihoodMap.size() == 0 )
                return null;

            for ( PerReadAlleleLikelihoodMap perReadLikelihoods : perReadAlleleLikelihoodMap.values() )
                totalSize += perReadLikelihoods.size();

            qualities = new int[totalSize];
            for ( PerReadAlleleLikelihoodMap perReadLikelihoods : perReadAlleleLikelihoodMap.values() ) {
                for (PileupElement p : perReadLikelihoods.getStoredPileupElements())
                    index = fillMappingQualitiesFromPileupAndUpdateIndex(p, index, qualities);


        }
        }
        else
            return null;



        double rms = MathUtils.rms(qualities);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", rms));
        return map;
    }

    private static int fillMappingQualitiesFromPileupAndUpdateIndex(final PileupElement p, final int inputIdx, final int[] qualities) {
        int outputIdx = inputIdx;
        if ( p.getMappingQual() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE )
            qualities[outputIdx++] = p.getMappingQual();

        return outputIdx;
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.RMS_MAPPING_QUALITY_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}