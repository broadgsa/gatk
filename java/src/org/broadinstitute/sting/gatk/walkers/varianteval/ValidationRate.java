package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
@Analysis(name = "Validation Rate", description = "Validation Rate")
public class ValidationRate extends VariantEvaluator {
    @DataPoint(name="# mono in comp",description = "Number of mono calls in the comparison ROD")
    long n_mono_in_comp;
    @DataPoint(name="# poly in comp",description = "Number of poly calls in the comparison ROD")
    long n_poly_in_comp;
    @DataPoint(name="% poly in comp",description = "Percent of poly calls in the comparison ROD")
    double percent_poly_in_comp;
    @DataPoint(name="# mono calls at mono sites",description = "Number of mono calls at mono sites")
    long n_mono_calls_at_mono_sites;
    @DataPoint(name="# poly calls at mono sites",description = "Number of poly calls at mono sites")
    long n_poly_calls_at_mono_sites;
    @DataPoint(name="# nocalls at mono sites",description = "Number of no calls at mono sites")
    long n_nocalls_at_mono_sites;
    @DataPoint(name="# mono sites called poly",description = "Percentage of mono sites called poly")
    double percent_mono_sites_called_poly;
    @DataPoint(name="# mono calls at poly sites",description = "Number of mono calls at poly sites")
    long n_mono_calls_at_poly_sites;
    @DataPoint(name="# poly calls at poly sites",description = "Number of poly calls at poly sites")
    long n_poly_calls_at_poly_sites;
    @DataPoint(name="# nocalls at poly sites",description = "Number of no calls at poly sites")
    long n_nocalls_at_poly_sites;
    @DataPoint(name="% poly sites called poly",description = "Number of poly sites called poly")
    double percent_poly_sites_called_poly;
    @DataPoint(description = "The PPV")
    double PPV;
    @DataPoint(description = "The sensitivity")
    double sensitivity;
    @DataPoint(description = "The False Discovery Rate (FDR)")
    double falseDiscoveryRate;

    // todo -- subset validation data by list of samples, if provided
    class SiteStats {
        long nPoly = 0;
        long nMono = 0;
        long nNoCall = 0;

        double polyPercent() {
            return 100 * rate(nPoly, nPoly + nMono + nNoCall);
        }
    }

    private SiteStats validationStats = new SiteStats();
    private SiteStats evalOverlapAtMono = new SiteStats();
    private SiteStats evalOverlapAtPoly = new SiteStats();

    public ValidationRate(VariantEvalWalker parent) {
        super(parent);
    }

    public String getName() {
        return "validationRate";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    @Override
    public void finalizeEvaluation() {
        long TP = evalOverlapAtPoly.nPoly; //  + evalOverlapAtMono.nMono + evalOverlapAtMono.nNoCall;
        long FP = evalOverlapAtMono.nPoly; //  + evalOverlapAtPoly.nMono;
        long FN = evalOverlapAtPoly.nMono + evalOverlapAtPoly.nNoCall;

        // fill in the output fields
        n_mono_in_comp = validationStats.nMono;
        n_poly_in_comp = validationStats.nPoly;
        percent_poly_in_comp = validationStats.polyPercent();

        n_mono_calls_at_mono_sites = evalOverlapAtMono.nMono;
        n_poly_calls_at_mono_sites = evalOverlapAtMono.nPoly;
        n_nocalls_at_mono_sites = evalOverlapAtMono.nNoCall;
        percent_mono_sites_called_poly = evalOverlapAtMono.polyPercent();

        n_mono_calls_at_poly_sites = evalOverlapAtPoly.nMono;
        n_poly_calls_at_poly_sites = evalOverlapAtPoly.nPoly;
        n_nocalls_at_poly_sites = evalOverlapAtPoly.nNoCall;
        percent_poly_sites_called_poly = evalOverlapAtPoly.polyPercent();
        PPV = 100 * rate(TP, TP + FP);
        sensitivity = 100 * rate(TP, TP + FN);
        falseDiscoveryRate = 100 * rate(FP, FP + TP);
    }

    public boolean enabled() {
        return true;
    }

    public String update2(VariantContext eval, VariantContext rawValidationData, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String interesting = null;

        if (rawValidationData!= null && rawValidationData.hasGenotypes() && rawValidationData.isNotFiltered()) {
            VariantContext validation = rawValidationData;

            SiteStats overlap;

            if (validation.isPolymorphic()) {
                validationStats.nPoly++;
                overlap = evalOverlapAtPoly;
                if (eval == null || eval.isMonomorphic())
                    interesting = "ValidationStatus=FN";
            } else {
                validationStats.nMono++;
                overlap = evalOverlapAtMono;

                if (eval != null && eval.isPolymorphic())
                    interesting = "ValidationStatus=FP";
            }

            if (eval == null)
                overlap.nNoCall++;
            //else if (eval.isPolymorphic()) // requires genotypes
            else if (eval.isSNP())
                overlap.nPoly++;
            else
                overlap.nMono++;
        }

        return interesting; // we don't capture any interesting sites
    }
}