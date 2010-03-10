package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class ValidationRate extends VariantEvaluator {
    class SiteStats {
        long nPoly = 0, nMono = 0, nNoCall = 0;

        double polyPercent() { return rate(nPoly, nPoly + nMono + nNoCall); }
    }

    private SiteStats validationStats = new SiteStats();
    private SiteStats evalOverlapAtMono = new SiteStats();
    private SiteStats evalOverlapAtPoly = new SiteStats();

    public ValidationRate(VariantEval2Walker parent) {
        // don't do anything
    }

    public String getName() {
        return "validationRate";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public String toString() {
        return getName() + ": " + summaryLine();
    }

    private String summaryLine() {
        return String.format("%d %d %.2f %d %d %d %.2f %d %d %d %.2f",
                validationStats.nMono, validationStats.nPoly, validationStats.polyPercent(),
                evalOverlapAtMono.nMono, evalOverlapAtMono.nPoly, evalOverlapAtMono.nNoCall, evalOverlapAtMono.polyPercent(),
                evalOverlapAtPoly.nMono, evalOverlapAtPoly.nPoly, evalOverlapAtPoly.nNoCall, evalOverlapAtPoly.polyPercent());
    }

    private static List<String> HEADER =
            Arrays.asList("n_mono_in_comp", "n_poly_in_comp", "percent_poly_in_comp",
                    "n_mono_calls_at_mono_sites", "n_poly_calls_at_mono_sites", "n_nocalls_at_mono_sites", "percent_mono_sites_called_poly",
                    "n_mono_calls_at_poly_sites", "n_poly_calls_at_poly_sites", "n_nocalls_at_poly_sites", "percent_poly_sites_called_poly");

    // making it a table
    public List<String> getTableHeader() {
        return HEADER;
    }

    public boolean enabled() { return true; }

    public List<List<String>> getTableRows() {
        return Arrays.asList(Arrays.asList(summaryLine().split(" ")));
    }

    public String update2(VariantContext eval, VariantContext validation, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String interesting = null;

        if ( validation != null && validation.hasGenotypes() && validation.isNotFiltered() ) {
            SiteStats overlap = null;

            if ( validation.isPolymorphic() ) {
                validationStats.nPoly++;
                overlap = evalOverlapAtPoly;
                if ( eval == null || eval.isMonomorphic() )
                    interesting = "ValidationStatus=FN";
            }
            else {
                validationStats.nMono++;
                overlap = evalOverlapAtMono;

                if ( eval != null && eval.isPolymorphic() )
                    interesting = "ValidationStatus=FP";
            }

            if ( eval == null )
                overlap.nNoCall++;
            else if ( eval.isPolymorphic() )
                overlap.nPoly++;
            else
                overlap.nMono++;
        }

        return interesting; // we don't capture any interesting sites
    }
}