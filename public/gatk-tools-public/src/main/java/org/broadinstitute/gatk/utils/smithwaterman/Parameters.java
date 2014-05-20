/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.gatk.utils.smithwaterman;

/**
 * Holds the core Smith-Waterman alignment parameters of
 *
 * match value, and mismatch, gap open and gap extension penalties
 *
 * User: depristo
 * Date: 4/11/13
 * Time: 12:03 PM
 */
public final class Parameters {
    public final double w_match;
    public final double w_mismatch;
    public final double w_open;
    public final double w_extend;
    public final double epsilon;

    /**
     * Create a new set of SW parameters
     * @param w_match the match score
     * @param w_mismatch the mismatch penalty
     * @param w_open the gap open penalty
     * @param w_extend the gap extension penalty
     * @param epsilon weight comparison error
     */
    public Parameters(final double w_match, final double w_mismatch, final double w_open, final double w_extend, final double epsilon) {
        if ( w_mismatch > 0 ) throw new IllegalArgumentException("w_mismatch must be <= 0 but got " + w_mismatch);
        if ( w_open> 0 ) throw new IllegalArgumentException("w_open must be <= 0 but got " + w_open);
        if ( w_extend> 0 ) throw new IllegalArgumentException("w_extend must be <= 0 but got " + w_extend);
        if ( Double.isNaN(epsilon)) throw new IllegalArgumentException("epsilon cannot be a NaN");
        if ( epsilon < 0.0 ) throw new IllegalArgumentException("epsilon cannot be negative");

        this.w_match = w_match;
        this.w_mismatch = w_mismatch;
        this.w_open = w_open;
        this.w_extend = w_extend;
        this.epsilon = epsilon;
    }

}
