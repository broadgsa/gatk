package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;


public class VECIndelArtifact implements VariantExclusionCriterion {
    private boolean exclude;
    private String source = "N/A";

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
        }
    }

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
        // TODO - fix indel call capability to span full indel
        if ( indelCall != null ) {
            exclude = true;
            source = "IndelCall";
            return;
        }

        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        // TODO - fix dbsnp capability to span full indel
        if ( dbsnp != null && dbsnp.isIndel() ) {
            exclude = true;
            source = "dbsnp";
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

    public boolean useZeroQualityReads() { return false; }
}