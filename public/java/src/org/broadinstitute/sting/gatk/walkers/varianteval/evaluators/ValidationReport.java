package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Collection;
import java.util.Set;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
@Analysis(description = "Assess site accuracy and sensitivity of callset against follow-up validation assay")
public class ValidationReport extends VariantEvaluator implements StandardEval {
    // todo -- note this isn't strictly allele away.  It's really focused on sites.  A/T call at a validated A/G site is currently counted as a TP
    @DataPoint(description = "nComp") int nComp = 0;
    @DataPoint(description = "TP") int TP = 0;
    @DataPoint(description = "FP") int FP = 0;
    @DataPoint(description = "FN") int FN = 0;
    @DataPoint(description = "TN") int TN = 0;

    @DataPoint(description = "Sensitivity") double sensitivity = 0;
    @DataPoint(description = "Specificity") double specificity = 0;
    @DataPoint(description = "PPV") double PPV = 0;
    @DataPoint(description = "FDR") double FDR = 0;

    @DataPoint(description = "CompMonoEvalNoCall") int CompMonoEvalNoCall = 0;
    @DataPoint(description = "CompMonoEvalFiltered") int CompMonoEvalFiltered = 0;
    @DataPoint(description = "CompMonoEvalMono") int CompMonoEvalMono = 0;
    @DataPoint(description = "CompMonoEvalPoly") int CompMonoEvalPoly = 0;

    @DataPoint(description = "CompPolyEvalNoCall") int CompPolyEvalNoCall = 0;
    @DataPoint(description = "CompPolyEvalFiltered") int CompPolyEvalFiltered = 0;
    @DataPoint(description = "CompPolyEvalMono") int CompPolyEvalMono = 0;
    @DataPoint(description = "CompPolyEvalPoly") int CompPolyEvalPoly = 0;

    @DataPoint(description = "CompFiltered") int CompFiltered = 0;
    @DataPoint(description = "Eval and comp have different alleles") int nDifferentAlleleSites = 0;

    private static final boolean TREAT_ALL_SITES_IN_EVAL_VCF_AS_CALLED = true;
    private static final boolean REQUIRE_IDENTICAL_ALLELES = false;

    private enum SiteStatus { NO_CALL, FILTERED, MONO, POLY }

    // Counts of ValidationSiteStatus x CallSiteStatus
    final int[][] counts = new int[SiteStatus.values().length][SiteStatus.values().length];

    @Override public int getComparisonOrder() { return 2; }
    @Override public boolean enabled() { return true; }

    @Override
    public void finalizeEvaluation() {
        for ( SiteStatus x : SiteStatus.values() )
            CompFiltered += getCounts(SiteStatus.FILTERED, x);

        CompMonoEvalNoCall = getCounts(SiteStatus.MONO, SiteStatus.NO_CALL);
        CompMonoEvalFiltered = getCounts(SiteStatus.MONO, SiteStatus.FILTERED);
        CompMonoEvalMono = getCounts(SiteStatus.MONO, SiteStatus.MONO);
        CompMonoEvalPoly = getCounts(SiteStatus.MONO, SiteStatus.POLY);

        CompPolyEvalNoCall = getCounts(SiteStatus.POLY, SiteStatus.NO_CALL);
        CompPolyEvalFiltered = getCounts(SiteStatus.POLY, SiteStatus.FILTERED);
        CompPolyEvalMono = getCounts(SiteStatus.POLY, SiteStatus.MONO);
        CompPolyEvalPoly = getCounts(SiteStatus.POLY, SiteStatus.POLY);

        TP = CompPolyEvalPoly;
        FN = CompPolyEvalNoCall + CompPolyEvalFiltered + CompPolyEvalMono;
        FP = CompMonoEvalPoly;
        TN = CompMonoEvalNoCall + CompMonoEvalFiltered + CompMonoEvalMono;

        for ( SiteStatus x : SiteStatus.values() )
            for ( SiteStatus y : SiteStatus.values() )
                nComp += getCounts(x, y);

        if ( nComp != TP + FN + FP + TN + CompFiltered )
            throw new ReviewedStingException("BUG: nComp != TP + FN + FP + TN + CompFiltered!");

        sensitivity = (100.0 * TP) / (TP + FN);
        specificity = (TN+FP > 0) ? (100.0 * TN) / (TN + FP) : 100.0;
        PPV = (100.0 * TP) / (TP + FP);
        FDR = (100.0 * FP) / (FP + TP);
    }

    private int getCounts(SiteStatus comp, SiteStatus eval) {
        return counts[comp.ordinal()][eval.ordinal()];
    }

    @Override
    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( comp != null ) { // we only need to consider sites in comp
            if ( REQUIRE_IDENTICAL_ALLELES && (eval != null && haveDifferentAltAlleles(eval, comp)))
                nDifferentAlleleSites++;
            else {
                SiteStatus evalStatus = calcSiteStatus(eval);
                SiteStatus compStatus = calcSiteStatus(comp);
                counts[compStatus.ordinal()][evalStatus.ordinal()]++;
            }
        }

        return null; // we don't capture any interesting sites
    }

    //
    // helper routines
    //
    public SiteStatus calcSiteStatus(VariantContext vc) {
        if ( vc == null ) return SiteStatus.NO_CALL;
        if ( vc.isFiltered() ) return SiteStatus.FILTERED;
        if ( vc.isMonomorphic() ) return SiteStatus.MONO;
        if ( vc.hasGenotypes() ) return SiteStatus.POLY;  // must be polymorphic if isMonomorphic was false and there are genotypes

        if ( vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
            int ac = 0;
            if ( vc.getNAlleles() > 2 ) {
                return SiteStatus.POLY;
////                System.out.printf("multiple alleles %s = %s%n", vc.getAlleles(), vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
//                // todo -- omg this is painful.  We need a better approach to dealing with multi-valued attributes
//                for ( String v : (List<String>)vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY) )
//                    ac += Integer.valueOf(v);
////                System.out.printf("  ac = %d%n", ac);
            }
            else
                ac = vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0);
            return ac > 0 ? SiteStatus.POLY : SiteStatus.MONO;
        } else {
            return TREAT_ALL_SITES_IN_EVAL_VCF_AS_CALLED ? SiteStatus.POLY : SiteStatus.NO_CALL; // we can't figure out what to do
            //return SiteStatus.NO_CALL; // we can't figure out what to do
        }
    }



    public boolean haveDifferentAltAlleles(VariantContext eval, VariantContext comp) {
        Collection<Allele> evalAlts = eval.getAlternateAlleles();
        Collection<Allele> compAlts = comp.getAlternateAlleles();
        if ( evalAlts.size() != compAlts.size() ) {
            return true;
        } else {
            // same size => every alt from eval must be in comp
            for ( Allele a : evalAlts ) {
                if ( ! compAlts.contains(a) ) {
//                    System.out.printf("Different alleles: %s:%d eval=%s comp=%s\n\t\teval=%s\n\t\tcomp=%s%n",
//                            eval.getChr(), eval.getStart(), eval.getAlleles(), comp.getAlleles(), eval, comp);
                    return true;
                }
            }

            return false;
        }
    }
}
