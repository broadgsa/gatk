package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.utils.MathUtils;
import net.sf.samtools.SAMRecord;

import java.util.List;


public class VECMappingQualityZero implements VariantExclusionCriterion {
    private int maximum = 50;
    private int mq0Count;
    private boolean exclude;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            maximum = Integer.valueOf(arguments);
        }
    }

    public void compute(char ref, AlignmentContext context, rodVariants variant) {
        List<SAMRecord> reads = context.getReads();
        mq0Count = 0;
        for (int i=0; i < reads.size(); i++) {
            if ( reads.get(i).getMappingQuality() == 0 )
                mq0Count++;
        }
        exclude = mq0Count > maximum;
    }

    public boolean isExcludable() {
        return exclude;
    }

    public String getStudyHeader() {
        return "MappingQualityZero("+maximum+")\tMQ0_count";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + mq0Count;
    }

    public boolean useZeroQualityReads() { return true; }
}