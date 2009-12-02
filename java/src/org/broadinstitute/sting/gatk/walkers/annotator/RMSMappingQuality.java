package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import net.sf.samtools.SAMRecord;

import java.util.List;


public class RMSMappingQuality extends StandardVariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {
        List<SAMRecord> reads = pileup.getReads();
        int[] qualities = new int[reads.size()];
        for (int i=0; i < reads.size(); i++)
            qualities[i] = reads.get(i).getMappingQuality();
        double rms = MathUtils.rms(qualities);
        return new Pair<String, String>(getKeyName(), String.format("%.2f", rms));
    }

    public String getKeyName() { return "MQ"; }

    public String getDescription() { return "MQ,1,Float,\"RMS Mapping Quality\""; }

    public boolean useZeroQualityReads() { return true; }
}