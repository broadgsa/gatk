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

package org.broadinstitute.gatk.engine.filters;

import htsjdk.samtools.SAMRecord;

/**
 * Filter out reads that are not paired, have their mate unmapped, are duplicates, fail vendor quality check or both mate and read are in the same strand.
 *
 * @author chartl
 * @since 5/18/11
 */
public class MateSameStrandFilter extends ReadFilter {

    public boolean filterOut(SAMRecord read) {
        return (! read.getReadPairedFlag() ) || read.getMateUnmappedFlag() || read.getDuplicateReadFlag() ||
                read.getReadFailsVendorQualityCheckFlag() || (read.getMateNegativeStrandFlag() == read.getReadNegativeStrandFlag());
    }
}
