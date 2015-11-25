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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Tag this stratification as dynamically determining the final strat based on the input data
 *
 * The paradigm here is simple.  We upfront create a strat with N states that reflect the finest grained
 * possible division of the data.  The data is processed, and statistics collected for each of the N states.
 * An update call is made to the stratification for evaluation VariantContext during each map call,
 * allowing the strat to collect data about the usage of each state.  A final call requests that
 * the stratification map down the N states into M states (typically less than N, not necessarily
 * a subset of N).  This is provided by returning a map from each of M state -> N states and
 * the VariantEval walker will combine all of the evaluations for N into a single value for
 * each M.
 *
 * For example, suppose I have a dynamic strat called AC, adopting 7 possible values 0,1,2,3,4,5,6.  This
 * strats tracks the number of eval vcs for each state, with final counts 0=1, 1=100, 2=10, 3=5, 4=3, 5=2, 6=1.
 * The stratification attempts to combine the strats down to so that each state has approximately the same
 * fraction of the data in each bin.  Overall there is 1+100+10+5+3+2+1=124 observations and 7 bins so we really
 * want ~ 18 observations in each bin.  So we merge 3-6 with 5+3+2+1 = 11 and keep 2, 1, and 0 as distinct bins.  We
 * return a map from 0 -> 0, 1 -> 1, 2 -> 2, 3-6 -> {3,4,5,6}.
 *
 * TODO - some open implementation questions
 * -- We should only create one stratifier overall.  How do we track this?  When we create the stratifiers
 *    perhaps we can look at them and create a tracker?
 * -- How do we create a new stratifier based on the finalStratifications() given the framework?  Conceptually
 *    this new thing is itself a stratifier, just like before, but it's states are determined at the end.  We'd
 *    then like to call not getRelevantStates but a different function that accepts an old state and returns
 *    the new state.  Perhaps the process should look like:
 *          finalizeStratification -> new Stratifier whose states are the final ones
 *          getNewState(old state) -> new state (one of those in getFinalStratification)
 *
 * @author Mark DePristo
 * @since 4/9/12
 */
public interface DynamicStratification {
    public void update(final VariantContext eval);
    public VariantStratifier finalizeStratification();
    public Object getFinalState(final Object oldState);
}
