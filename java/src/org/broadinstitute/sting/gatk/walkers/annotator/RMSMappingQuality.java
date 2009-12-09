package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;


public class RMSMappingQuality extends StandardVariantAnnotation {

    public String annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation) {
        int[] qualities = new int[pileup.size()];
        int index = 0;
        for (PileupElement p : pileup )
            qualities[index++] = p.getRead().getMappingQuality();
        double rms = MathUtils.rms(qualities);
        return String.format("%.2f", rms);
    }

    public String getKeyName() { return VCFRecord.RMS_MAPPING_QUALITY_KEY; }

    public String getDescription() { return getKeyName() + ",1,Float,\"RMS Mapping Quality\""; }

    public boolean useZeroQualityReads() { return true; }
}