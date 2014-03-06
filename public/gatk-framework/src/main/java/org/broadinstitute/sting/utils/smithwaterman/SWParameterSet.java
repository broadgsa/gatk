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

package org.broadinstitute.sting.utils.smithwaterman;

/**
 * Handy named collection of common Smith-waterman parameter sets
 *
 * User: depristo
 * Date: 4/11/13
 * Time: 12:02 PM
 */
public enum SWParameterSet {
    // match=1, mismatch = -1/3, gap=-(1+k/3)
    ORIGINAL_DEFAULT(new Parameters(1.0,-1.0/3.0,-1.0-1.0/3.0,-1.0/3.0)),

    /**
     * A standard set of values for NGS alignments
     */
    STANDARD_NGS(new Parameters(5.0, -10.0, -22.0, -1.2));

    protected Parameters parameters;

    SWParameterSet(final Parameters parameters) {
        if ( parameters == null ) throw new IllegalArgumentException("parameters cannot be null");

        this.parameters = parameters;
    }
}
