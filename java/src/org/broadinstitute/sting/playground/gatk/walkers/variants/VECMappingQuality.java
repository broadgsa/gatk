package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import net.sf.samtools.SAMRecord;

import java.util.List;


public class VECMappingQuality implements VariantExclusionCriterion {
    private double minQuality = 5.0;
    private double rms;
    private boolean exclude;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            minQuality = Double.valueOf(arguments);
        }
    }

    public void compute(char ref, LocusContext context, rodVariants variant) {
        List<SAMRecord> reads = context.getReads();

        rms = 0.0;
        for (int readIndex = 0; readIndex < reads.size(); readIndex++) {
            int qual = reads.get(readIndex).getMappingQuality();
            rms += qual * qual;
        }
        rms /= reads.size();
        rms = Math.sqrt(rms);
    }

    public boolean isExcludable() {
        return exclude;
    }

    public String getStudyHeader() {
        return "MappingQuality("+minQuality+")\trms";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + rms;
    }

    public boolean useZeroQualityReads() { return true; }
}