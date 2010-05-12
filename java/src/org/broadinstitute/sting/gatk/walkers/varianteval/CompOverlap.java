package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
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
@Analysis(name = "Comp Overlap", description = "the overlap between eval and comp sites")
public class CompOverlap extends VariantEvaluator {

    @DataPoint(name = "eval sites", description = "number of eval SNP sites")
    long nEvalSNPs = 0;

    @DataPoint(name = "comp sites", description = "number of comp SNP sites")
    long nCompSNPs = 0;

    @DataPoint(name = "evals not at comp", description = "number of eval sites outside of comp sites")
    long novelSites = 0;

    @DataPoint(name = "evals at comp", description = "number of eval sites at comp sites")
    long nSNPsAtComp = 0;

    @DataPoint(name = "% evals at comp", description = "percentage of eval sites at comp sites")
    double compRate = 0.0;

    @DataPoint(name = "concordant", description = "number of concordant sites")
    long nConcordant = 0;

    @DataPoint(name = "% concordant", description = "the concordance rate")
    double concordantRate = 0.0;

    public CompOverlap(VariantEvalWalker parent) {
        super(parent);
    }

    public String getName() {
        return "compOverlap";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public long nNovelSites() { return nEvalSNPs - nSNPsAtComp; }
    public double compRate() { return rate(nSNPsAtComp, nEvalSNPs); }
    public double concordanceRate() { return rate(nConcordant, nSNPsAtComp); }

    public void finalizeEvaluation() {
        compRate = 100 * compRate();
        concordantRate = 100 * concordanceRate();
        novelSites = nNovelSites();
    }

    public boolean enabled() {
        return true;
    }

    /**
     * Returns true if every allele in eval is also in comp
     *
     * @param eval  eval context
     * @param comp db context
     * @return true if eval and db are discordant
     */
    public boolean discordantP(VariantContext eval, VariantContext comp) {
        for (Allele a : eval.getAlleles()) {
            if (!comp.hasAllele(a, true))
                return true;
        }

        return false;
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        boolean compIsGood = comp != null && comp.isSNP() && comp.isNotFiltered();
        boolean evalIsGood = eval != null && eval.isSNP();

        if (compIsGood) nCompSNPs++;           // count the number of comp events
        if (evalIsGood) nEvalSNPs++;           // count the number of eval events

        if (compIsGood && evalIsGood) {
            nSNPsAtComp++;

            if (!discordantP(eval, comp))     // count whether we're concordant or not with the comp value
                nConcordant++;
        }

        return null; // we don't capture any interesting sites
    }
}
