/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * Created by IntelliJ IDEA.
 * User: depristo, ebanks
 * Date: Jan 22, 2010
 * Time: 2:25:19 PM
 *
 * Useful helper class to communicate the results of calculateGenotype to framework
 */
public class VariantCallContext extends VariantContext {

    // Was the site called confidently, either reference or variant?
    public boolean confidentlyCalled = false;

    // Should this site be emitted?
    public boolean shouldEmit = true;

    VariantCallContext(VariantContext vc, boolean confidentlyCalledP) {
        super(vc);
        this.confidentlyCalled = confidentlyCalledP;
    }

    VariantCallContext(VariantContext vc, boolean confidentlyCalledP, boolean shouldEmit) {
        super(vc);
        this.confidentlyCalled = confidentlyCalledP;
        this.shouldEmit = shouldEmit;
    }

    /* these methods are only implemented for GENOTYPE_GIVEN_ALLELES MODE */
    //todo -- expand these methods to all modes

    /**
     *
     * @param callConfidenceThreshold the Unified Argument Collection STANDARD_CONFIDENCE_FOR_CALLING
     * @return true if call was confidently ref
     */
    public boolean isCalledRef(double callConfidenceThreshold) {
        return (confidentlyCalled && (getPhredScaledQual() < callConfidenceThreshold));
    }

    /**
     *
     * @param callConfidenceThreshold the Unified Argument Collection STANDARD_CONFIDENCE_FOR_CALLING
     * @return true if call was confidently alt
     */
    public boolean isCalledAlt(double callConfidenceThreshold) {
        return (confidentlyCalled && (getPhredScaledQual() > callConfidenceThreshold));
    }

}