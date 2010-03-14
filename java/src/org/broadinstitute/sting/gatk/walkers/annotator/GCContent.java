package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;


public class GCContent implements VariantAnnotation {

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        double content = computeGCContent(ref);
        return String.format("%.2f", content);
    }

    public String getKeyName() { return "GC"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("GC", 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "GC content within 20 bp +/- the variant"); }

    public boolean useZeroQualityReads() { return false; }

    private static double computeGCContent(ReferenceContext ref) {
        int gc = 0, at = 0;

        for ( char base : ref.getBases() ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
            if ( baseIndex == BaseUtils.gIndex || baseIndex == BaseUtils.cIndex )
                gc++;
            else if ( baseIndex == BaseUtils.aIndex || baseIndex == BaseUtils.tIndex )
                at++;
            else
                ; // ignore
        }

        int sum = gc + at;
        return (100.0*gc) / (sum == 0 ? 1 : sum);
     }
}