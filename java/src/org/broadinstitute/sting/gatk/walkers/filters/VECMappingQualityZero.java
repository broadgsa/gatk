package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.VariantContext;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.HashMap;


public class VECMappingQualityZero implements VariantExclusionCriterion {
    private int maximum = 50;
    private int mq0Count;
    private boolean exclude;

    public void initialize(HashMap<String,String> args) {
        if ( args.get("max") != null )
            maximum = Integer.valueOf(args.get("max"));
    }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        List<SAMRecord> reads = context.getAlignmentContext(useZeroQualityReads()).getReads();
        mq0Count = 0;
        for (int i=0; i < reads.size(); i++) {
            if ( reads.get(i).getMappingQuality() == 0 )
                mq0Count++;
        }
        exclude = mq0Count > maximum;
    }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return "MappingQualityZero("+maximum+")\tMQ0_count";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + mq0Count;
    }

    public String getVCFFilterString() {
        return "MQzero" + mq0Count;
    }

    public boolean useZeroQualityReads() { return true; }
}