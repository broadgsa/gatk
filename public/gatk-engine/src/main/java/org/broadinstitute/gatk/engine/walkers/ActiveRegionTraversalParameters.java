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

package org.broadinstitute.gatk.engine.walkers;

import org.broadinstitute.gatk.utils.activeregion.BandPassActivityProfile;

import java.lang.annotation.Documented;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Describes the parameters that this walker requires of the active region traversal
 *
 * User: rpoplin
 * Date: 1/18/12
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)

public @interface ActiveRegionTraversalParameters {
    /**
     * How far to either side of the active region itself should we include reads?
     *
     * That is, if the active region is 10 bp wide, and extension is 5, ART will provide
     * the walker with active regions 10 bp, with 5 bp of extension on either side, and
     * all reads that cover the 20 bp of the region + extension.
     *
     * @return the size of the active region extension we'd like
     */
    public int extension() default 0;

    /**
     * The minimum number of bp for an active region, when we need to chop it up into pieces because
     * it's become too big.  This only comes into effect when there's literally no good place to chop
     * that does make the region smaller than this value.
     *
     * @return the min size in bp of regions
     */
    public int minRegion() default 50;

    /**
     * The maximum size in bp of active regions wanted by this walker
     *
     * Active regions larger than this value are automatically cut up by ART into smaller
     * regions of size <= this value.
     *
     * @return the max size in bp of regions
     */
    public int maxRegion() default 1500;

    /**
     * The variance value for the Gaussian kernel of the band pass filter employed by ART
     * @return the breadth of the band pass gaussian kernel we want for our traversal
     */
    public double bandPassSigma() default BandPassActivityProfile.DEFAULT_SIGMA;

    /**
     * What is the maximum number of reads we're willing to hold in memory per sample
     * during the traversal?  This limits our exposure to unusually large amounts
     * of coverage in the engine.
     * @return the maximum number of reads we're willing to hold in memory
     */
    public int maxReadsToHoldInMemoryPerSample() default 30000;

    /**
     * No matter what the per sample value says, we will never hold more than this
     * number of reads in memory at any time.  Provides an upper bound on the total number
     * of reads in the case where we have a lot of samples.
     * @return the maximum number of reads to hold in memory
     */
    public int maxReadsToHoldTotal() default 10000000;
}
