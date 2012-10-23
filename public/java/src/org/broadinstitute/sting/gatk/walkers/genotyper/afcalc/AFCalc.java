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

package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.util.List;


/**
 * Generic interface for calculating the probability of alleles segregating given priors and genotype likelihoods
 */
public abstract class AFCalc implements Cloneable {
    private final static Logger defaultLogger = Logger.getLogger(AFCalc.class);

    protected final int nSamples;
    protected final int maxAlternateAllelesToGenotype;

    protected Logger logger = defaultLogger;

    private SimpleTimer callTimer = new SimpleTimer();
    private final StateTracker stateTracker;
    private ExactCallLogger exactCallLogger = null;

    /**
     * Create a new AFCalc object capable of calculating the prob. that alleles are
     * segregating among nSamples with up to maxAltAlleles for SNPs and maxAltAllelesForIndels
     * for indels for samples with ploidy
     *
     * @param nSamples number of samples, must be > 0
     * @param maxAltAlleles maxAltAlleles for SNPs
     * @param ploidy the ploidy, must be > 0
     */
    protected AFCalc(final int nSamples, final int maxAltAlleles, final int ploidy) {
        if ( nSamples < 0 ) throw new IllegalArgumentException("nSamples must be greater than zero " + nSamples);
        if ( maxAltAlleles < 1 ) throw new IllegalArgumentException("maxAltAlleles must be greater than zero " + maxAltAlleles);
        if ( ploidy < 1 ) throw new IllegalArgumentException("ploidy must be > 0 but got " + ploidy);

        this.nSamples = nSamples;
        this.maxAlternateAllelesToGenotype = maxAltAlleles;
        this.stateTracker = new StateTracker(maxAltAlleles);
    }

    /**
     * Enable exact call logging to file
     *
     * @param exactCallsLog the destination file
     */
    public void enableProcessLog(final File exactCallsLog) {
        exactCallLogger = new ExactCallLogger(exactCallsLog);
    }

    /**
     * Use this logger instead of the default logger
     *
     * @param logger
     */
    public void setLogger(Logger logger) {
        this.logger = logger;
    }

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information
     * @param log10AlleleFrequencyPriors a prior vector nSamples x 2 in length indicating the Pr(AF = i)
     * @return result (for programming convenience)
     */
    public AFCalcResult getLog10PNonRef(final VariantContext vc, final double[] log10AlleleFrequencyPriors) {
        if ( vc == null ) throw new IllegalArgumentException("VariantContext cannot be null");
        if ( log10AlleleFrequencyPriors == null ) throw new IllegalArgumentException("priors vector cannot be null");
        if ( stateTracker == null ) throw new IllegalArgumentException("Results object cannot be null");

        // reset the result, so we can store our new result there
        stateTracker.reset();

        final VariantContext vcWorking = reduceScope(vc);

        callTimer.start();
        final AFCalcResult result = computeLog10PNonRef(vcWorking, log10AlleleFrequencyPriors);
        final long nanoTime = callTimer.getElapsedTimeNano();

        if ( exactCallLogger != null )
            exactCallLogger.printCallInfo(vcWorking, log10AlleleFrequencyPriors, nanoTime, result);

        return result;
    }

    /**
     * Convert the final state of the state tracker into our result as an AFCalcResult
     *
     * Assumes that stateTracker has been updated accordingly
     *
     * @param vcWorking the VariantContext we actually used as input to the calc model (after reduction)
     * @param log10AlleleFrequencyPriors the priors by AC vector
     * @return a AFCalcResult describing the result of this calculation
     */
    @Requires("stateTracker.getnEvaluations() >= 0")
    @Ensures("result != null")
    protected AFCalcResult getResultFromFinalState(final VariantContext vcWorking, final double[] log10AlleleFrequencyPriors) {
        stateTracker.setAllelesUsedInGenotyping(vcWorking.getAlleles());
        return stateTracker.toAFCalcResult(log10AlleleFrequencyPriors);
    }

    // ---------------------------------------------------------------------------
    //
    // Abstract methods that should be implemented by concrete implementations
    // to actually calculate the AF
    //
    // ---------------------------------------------------------------------------

    /**
     * Look at VC and perhaps return a new one of reduced complexity, if that's necessary
     *
     * Used before the call to computeLog10PNonRef to simply the calculation job at hand,
     * if vc exceeds bounds.  For example, if VC has 100 alt alleles this function
     * may decide to only genotype the best 2 of them.
     *
     * @param vc the initial VC provided by the caller to this AFcalculation
     * @return a potentially simpler VC that's more tractable to genotype
     */
    @Requires("vc != null")
    @Ensures("result != null")
    protected abstract VariantContext reduceScope(final VariantContext vc);

    /**
     * Actually carry out the log10PNonRef calculation on vc, storing results in results
     *
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param log10AlleleFrequencyPriors        priors
     * @return a AFCalcResult object describing the results of this calculation
     */
    @Requires({"vc != null", "log10AlleleFrequencyPriors != null"})
    protected abstract AFCalcResult computeLog10PNonRef(final VariantContext vc,
                                                        final double[] log10AlleleFrequencyPriors);

    /**
     * Subset VC to the just allelesToUse, updating genotype likelihoods
     *
     * Must be overridden by concrete subclasses
     *
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes
     * @param ploidy
     * @return GenotypesContext object
     */
    public abstract GenotypesContext subsetAlleles(final VariantContext vc,
                                                   final List<Allele> allelesToUse,
                                                   final boolean assignGenotypes,
                                                   final int ploidy);

    // ---------------------------------------------------------------------------
    //
    // accessors
    //
    // ---------------------------------------------------------------------------

    public int getMaxAltAlleles() {
        return maxAlternateAllelesToGenotype;
    }

    protected StateTracker getStateTracker() {
        return stateTracker;
    }

}