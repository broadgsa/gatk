package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.utils.MathUtils;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.HashMap;


public class VECMappingQuality implements VariantExclusionCriterion {
    private double minQuality = 5.0;
    private double rms;
    private boolean exclude;

    public void initialize(HashMap<String,String> args) {
        if ( args.get("min") != null )
            minQuality = Double.valueOf(args.get("min"));
    }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        List<SAMRecord> reads = context.getAlignmentContext(useZeroQualityReads()).getReads();
        int[] qualities = new int[reads.size()];
        for (int i=0; i < reads.size(); i++)
            qualities[i] = reads.get(i).getMappingQuality();
        rms = MathUtils.rms(qualities);
        exclude = rms < minQuality;        
    }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return "MappingQuality("+minQuality+")\trms";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + rms;
    }

    public boolean useZeroQualityReads() { return true; }
}