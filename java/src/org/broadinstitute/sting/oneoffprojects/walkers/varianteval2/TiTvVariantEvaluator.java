package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.oneoffprojects.variantcontext.VariantContext;
import org.broadinstitute.sting.oneoffprojects.variantcontext.Genotype;
import org.broadinstitute.sting.utils.StingException;

import java.util.List;
import java.util.Arrays;

public class TiTvVariantEvaluator extends VariantEvaluator {
    long nTi = 0, nTv = 0;
    long nTiInStd = 0, nTvInStd = 0;

    public TiTvVariantEvaluator(VariantEval2Walker parent) {
        // don't do anything
    }

    public String getName() {
        return "titv";
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    public void updateTiTv(VariantContext vc, boolean updateStandard) {
        if ( vc != null && vc.isSNP() && vc.isBiallelic() ) {
            if ( vc.isTransition() ) {
                if ( updateStandard ) nTiInStd++; else nTi++;
            } else {
                if ( updateStandard ) nTvInStd++; else nTv++;
            }
        }
    }

    public void update2(VariantContext vc1, VariantContext vc2, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( vc1 != null ) updateTiTv(vc1, false);
        if ( vc2 != null ) updateTiTv(vc2, true);
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