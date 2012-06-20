package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
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

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Root Mean Square of the mapping quality of the reads across all samples.
 */
public class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int totalSize = 0;
        for ( AlignmentContext context : stratifiedContexts.values() )
            totalSize += context.size();

        final int[] qualities = new int[totalSize];
        int index = 0;

        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            AlignmentContext context = sample.getValue();
            final ReadBackedPileup pileup = context.getBasePileup();
            for (PileupElement p : pileup ) {
                if ( p.getMappingQual() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE )
                    qualities[index++] = p.getMappingQual();
            }
        }

        double rms = MathUtils.rms(qualities);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", rms));
        return map;
    }

    public Map<String, Object> annotate(Map<String, Map<Allele, List<GATKSAMRecord>>> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int depth = 0;
        for ( final Map<Allele, List<GATKSAMRecord>> alleleBins : stratifiedContexts.values() ) {
            for ( final Map.Entry<Allele, List<GATKSAMRecord>> alleleBin : alleleBins.entrySet() ) {
                depth += alleleBin.getValue().size();
            }
        }

        final int[] qualities = new int[depth];
        int index = 0;

        for ( final Map<Allele, List<GATKSAMRecord>> alleleBins : stratifiedContexts.values() ) {
            for ( final List<GATKSAMRecord> reads : alleleBins.values() ) {
                for ( final GATKSAMRecord read : reads ) {
                    if ( read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE )
                        qualities[index++] = read.getMappingQuality();
                }
            }
        }

        final Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", MathUtils.rms(qualities)));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFConstants.RMS_MAPPING_QUALITY_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}