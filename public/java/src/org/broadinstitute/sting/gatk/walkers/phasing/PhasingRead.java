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
package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Arrays;

public class PhasingRead extends BaseArray {
    private PreciseNonNegativeDouble mappingProb; // the probability that this read is mapped correctly
    private PreciseNonNegativeDouble[] baseProbs; // the probabilities that the base identities are CORRECT
    private PreciseNonNegativeDouble[] baseErrorProbs; // the probabilities that the base identities are INCORRECT

    public PhasingRead(int length, int mappingQual) {
        super(length);

        this.mappingProb = new PreciseNonNegativeDouble(QualityUtils.qualToProb((byte)mappingQual));

        this.baseProbs = new PreciseNonNegativeDouble[length];
        Arrays.fill(this.baseProbs, null);

        this.baseErrorProbs = new PreciseNonNegativeDouble[length];
        Arrays.fill(this.baseErrorProbs, null);
    }

    public void updateBaseAndQuality(int index, Byte base, byte baseQual) {
        updateBase(index, base);

        double errProb = QualityUtils.qualToErrorProb(baseQual);

        // The base error should be AT LEAST AS HIGH as the mapping error [equivalent to capping the base quality (BQ) by the mapping quality (MQ)]:
        errProb = Math.max(errProb, 1.0 - mappingProb.getValue());

        baseProbs[index] = new PreciseNonNegativeDouble(1.0 - errProb); // The probability that the true base is the base called in the read
        baseErrorProbs[index] = new PreciseNonNegativeDouble(errProb / 3.0); // DIVIDE up the error probability EQUALLY over the 3 non-called bases
    }

    public PhasingScore matchHaplotypeClassScore(HaplotypeClass hapClass) {
        PreciseNonNegativeDouble value = new PreciseNonNegativeDouble(0.0);
        for (Haplotype h : hapClass)
            value.plusEqual(matchHaplotypeScore(h));

        return new PhasingScore(value);
    }

    private PreciseNonNegativeDouble matchHaplotypeScore(Haplotype hap) {
        PreciseNonNegativeDouble score = new PreciseNonNegativeDouble(1.0);

        int sz = this.bases.length;
        if (sz != hap.bases.length)
            throw new ReviewedStingException("Read and Haplotype should have same length to be compared!");

        // Technically, this HAS NO EFFECT since it is multiplied in for ALL haplotype pairs, but do so for completeness:
        score.timesEqual(mappingProb);

        for (int i = 0; i < sz; i++) {
            Byte thisBase = this.getBase(i);
            Byte hapBase = hap.getBase(i);
            if (thisBase != null && hapBase != null) {
                if (BaseUtils.basesAreEqual(thisBase, hapBase))
                    score.timesEqual(baseProbs[i]);
                else
                    score.timesEqual(baseErrorProbs[i]);
            }
        }

        return score;
    }
}
