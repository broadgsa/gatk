/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import java.util.Set;
import java.util.HashSet;

/**
 * Class for maintaining the arguments that need to be passed to UnifiedGenotyper.runGenotyper()
 * (so that they only need to be computed one time)
 */
public class UGCalculationArguments {
    // should we annotate dbsnp?
    protected boolean annotateDbsnp = false;
    // should we annotate hapmap2?
    protected boolean annotateHapmap2 = false;
    // should we annotate hapmap3?
    protected boolean annotateHapmap3 = false;

    // the unified argument collection
    protected UnifiedArgumentCollection UAC = null;

    // the model used for calculating genotypes
    protected ThreadLocal<GenotypeCalculationModel> gcm = new ThreadLocal<GenotypeCalculationModel>();

    // samples in input
    protected Set<String> samples = new HashSet<String>();
}
