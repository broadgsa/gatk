/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.gatk.iterators.LocusIteratorFilter;

/**
 * Filters out read pairs where the reads are so long relative to the over fragment size that they are
 * reading into each other's adaptors.
 *
 * Normally, fragments are sufficiently far apart that reads aren't reading into each other.
 *
 * |-------------------->                                   first read
 *                                 <--------------------|   second read
 *
 * Sometimes, mostly due to lab errors or constraints, fragment library are made too short relative to the
 * length of the reads.  For example, it's possible to have 76bp PE reads with 125 bp inserts, so that ~25 bp of each
 * read overlaps with its mate.
 *
 * |-------------------->           first read
 *         <--------------------|   second read
 *
 * This filter deals with the situation where the fragment is so small that the each read actually reads into the
 * adaptor sequence of its mate, generating mismatches at both ends of the read:
 *
 *              |----------------XXXX>      first read
 *         <XXXX----------------|           second read
 *
 * So, how do we detect the problematic case?  Since the first and second pair can land on chromosome in either order, we need
 * to consider not first and second, but rather left and right, defined by which is :
 *
 *              |----------------XXXX>      left read
 *         <XXXX----------------|           right read
 *
 *
 * The problematic situation occurs when the alignment stop of this read is greater the alignment start of its mate
 *
 * @author depristo
 * @version 0.1
 */
public class InAdaptorFilter implements LocusIteratorFilter {
    public boolean filterOut(final SAMRecord rec, long basePos) {
        return ReadUtils.readPairBaseOverlapType(rec, basePos) == ReadUtils.OverlapType.IN_ADAPTOR;
    }
}