/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
@Analysis(description = "The overlap between eval and comp sites")
public class CompOverlap extends VariantEvaluator implements StandardEval {
    @DataPoint(description = "number of eval variant sites", format = "%d")
    public long nEvalVariants = 0;

    @DataPoint(description = "number of eval sites outside of comp sites", format = "%d")
    public long novelSites = 0;

    @DataPoint(description = "number of eval sites at comp sites", format = "%d")
    public long nVariantsAtComp = 0;

    @DataPoint(description = "percentage of eval sites at comp sites", format = "%.2f" )
    public double compRate = 0.0;

    @DataPoint(description = "number of concordant sites", format = "%d")
    public long nConcordant = 0;

    @DataPoint(description = "the concordance rate", format = "%.2f")
    public double concordantRate = 0.0;

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public long nNovelSites() { return nEvalVariants - nVariantsAtComp; }
    public double compRate() { return rate(nVariantsAtComp, nEvalVariants); }
    public double concordanceRate() { return rate(nConcordant, nVariantsAtComp); }

    public void finalizeEvaluation() {
        compRate = 100 * compRate();
        concordantRate = 100 * concordanceRate();
        novelSites = nNovelSites();
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

    public void update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        boolean evalIsGood = eval != null && eval.isPolymorphicInSamples();
        boolean compIsGood = comp != null && comp.isNotFiltered();

        if (evalIsGood) nEvalVariants++;           // count the number of eval events

        if (compIsGood && evalIsGood) {
            nVariantsAtComp++;

            if (!discordantP(eval, comp)) {    // count whether we're concordant or not with the comp value
                nConcordant++;
            }
        }
    }
}
