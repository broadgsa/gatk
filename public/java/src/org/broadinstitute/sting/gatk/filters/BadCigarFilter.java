/*
 * Copyright (c) 2009 The Broad Institute
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

package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import java.util.Iterator;

/**
 * Filter out reads with wonky cigar strings.
 *
 *  - No reads with Hard/Soft clips in the middle of the cigar
 *  - No reads starting with deletions (with or without preceding clips)
 *  - No reads ending in deletions (with or without follow-up clips)
 *  - No reads that are fully hard or soft clipped
 *  - No reads that have consecutive indels in the cigar (II, DD, ID or DI)
 *
 *  ps: apparently an empty cigar is okay...
 *
 * @author ebanks
 * @version 0.1
 */

public class BadCigarFilter extends ReadFilter {

    public boolean filterOut(final SAMRecord rec) {
        final Cigar c = rec.getCigar();

        // if there is no Cigar then it can't be bad
        if( c.isEmpty() ) {
            return false;
        }

        Iterator<CigarElement> elementIterator = c.getCigarElements().iterator();

        CigarOperator firstOp = CigarOperator.H;
        while (elementIterator.hasNext() && (firstOp == CigarOperator.H || firstOp == CigarOperator.S)) {
            CigarOperator op = elementIterator.next().getOperator();

            // No reads with Hard/Soft clips in the middle of the cigar
            if (firstOp != CigarOperator.H && op == CigarOperator.H) {
                    return true;
            }
            firstOp = op;
        }

        // No reads starting with deletions (with or without preceding clips)
        if (firstOp == CigarOperator.D) {
            return true;
        }

        boolean hasMeaningfulElements = (firstOp != CigarOperator.H && firstOp != CigarOperator.S);
        boolean previousElementWasIndel = firstOp == CigarOperator.I;
        CigarOperator lastOp = firstOp;
        CigarOperator previousOp = firstOp;

        while (elementIterator.hasNext()) {
            CigarOperator op = elementIterator.next().getOperator();

            if (op != CigarOperator.S && op != CigarOperator.H) {

                // No reads with Hard/Soft clips in the middle of the cigar
                if (previousOp == CigarOperator.S || previousOp == CigarOperator.H)
                    return true;

                lastOp = op;

                if (!hasMeaningfulElements && op.consumesReadBases()) {
                    hasMeaningfulElements = true;
                }

                if (op == CigarOperator.I || op == CigarOperator.D) {

                    // No reads that have consecutive indels in the cigar (II, DD, ID or DI)
                    if (previousElementWasIndel) {
                        return true;
                    }
                    previousElementWasIndel = true;
                }
                else {
                    previousElementWasIndel = false;
                }
            }
            // No reads with Hard/Soft clips in the middle of the cigar
            else if (op == CigarOperator.S && previousOp == CigarOperator.H) {
                return true;
            }

            previousOp = op;
        }

        // No reads ending in deletions (with or without follow-up clips)
        // No reads that are fully hard or soft clipped
        return lastOp == CigarOperator.D || !hasMeaningfulElements;
    }
}