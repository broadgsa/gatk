package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;

import java.util.List;
import java.util.Arrays;

@Analysis(name = "Ti/Tv Variant Evaluator", description = "Ti/Tv Variant Evaluator")
public class TiTvVariantEvaluator extends VariantEvaluator {

    @DataPoint(name = "ti_count", description = "number of transition loci")
    long nTi = 0;
    @DataPoint(name = "tv_count", description = "number of transversion loci")
    long nTv = 0;
    @DataPoint(name = "ti_count_std", description = "number of transition sites in the std")
    long nTiInStd = 0;
    @DataPoint(name = "tv_count_std", description = "number of transversion sites in the std")
    long nTvInStd = 0;

    public TiTvVariantEvaluator(VariantEval2Walker parent) {
        // don't do anything
    }

    public boolean enabled() {
        return true;
    }

    public String getName() {
        return "titv";
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    public void updateTiTv(VariantContext vc, boolean updateStandard) {
        if (vc != null && vc.isSNP() && vc.isBiallelic()) {
            if (vc.isTransition()) {
                if (updateStandard) nTiInStd++;
                else nTi++;
            } else {
                if (updateStandard) nTvInStd++;
                else nTv++;
            }
        }
    }

    public String update2(VariantContext vc1, VariantContext vc2, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (vc1 != null) updateTiTv(vc1, false);
        if (vc2 != null) updateTiTv(vc2, true);

        //if ( vc1 == null && vc2 != null && vc2.isSNP() && vc2.isBiallelic() )
        //    System.out.printf("VC2 = %s%n", vc2);
        //if ( vc2 != null && vc2.getName().equals("dbsnp") )

        return null; // we don't capture any intersting sites
    }

    public String toString() {
        return getName() + ": " + summaryLine();
    }

    private String summaryLine() {
        return String.format("%d %d %.2f %d %d %.2f",
                nTi, nTv, ratio(nTi, nTv),
                nTiInStd, nTvInStd, ratio(nTiInStd, nTvInStd));
    }

    private static List<String> HEADER =
            Arrays.asList("nTi", "nTv", "TiTvRatio", "nTiStandard", "nTvStandard", "TiTvRatioStandard");

    // making it a table
    public List<String> getTableHeader() {
        return HEADER;
    }

    public List<List<String>> getTableRows() {
        return Arrays.asList(Arrays.asList(summaryLine().split(" ")));
    }
}