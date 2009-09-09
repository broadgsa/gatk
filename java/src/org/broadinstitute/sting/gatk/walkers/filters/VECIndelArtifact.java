package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;

import java.util.HashMap;


public class VECIndelArtifact implements VariantExclusionCriterion {
    private boolean exclude;
    private String source = "N/A";

    public void initialize(HashMap<String,String> arguments) {}

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        RefMetaDataTracker tracker = context.getTracker();

        CleanedOutSNPROD cleanedSNP = (CleanedOutSNPROD)tracker.lookup("cleaned", null);
        if ( cleanedSNP != null && !cleanedSNP.isRealSNP() ) {
            exclude = true;
            source = "Cleaner";
            return;
        }

        AllelicVariant indelCall = (AllelicVariant)tracker.lookup("indels", null);
        if ( indelCall != null ) {
            exclude = true;
            source = "IndelCall";
            return;
        }

        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp != null && dbsnp.isIndel() ) {
            exclude = true;
            source = "dbSNP";
            return;
        }

        exclude = false;
    }

    public double inclusionProbability() {
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return "IndelArtifact\tSource";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + source;
    }

    public String getVCFFilterString() {
        return "InDel";
    }

    public boolean useZeroQualityReads() { return false; }
}