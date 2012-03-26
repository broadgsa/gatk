package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.BitSetUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.BitSet;

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

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 *
 * The Reported Quality Score covariate.
 */

public class QualityScoreCovariate implements RequiredCovariate {

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
    }

    @Override
    public CovariateValues getValues(final GATKSAMRecord read) {
        int readLength = read.getReadLength();

        BitSet[] mismatches = new BitSet[readLength];
        BitSet[] insertions = new BitSet[readLength];
        BitSet[] deletions = new BitSet[readLength];

        byte[] baseQualities = read.getBaseQualities();
        byte[] baseInsertionQualities = read.getBaseInsertionQualities();
        byte[] baseDeletionQualities = read.getBaseDeletionQualities();

        for (int i = 0; i < baseQualities.length; i++) {
            mismatches[i] = BitSetUtils.bitSetFrom(baseQualities[i]);
            insertions[i] = BitSetUtils.bitSetFrom(baseInsertionQualities[i]);
            deletions[i] = BitSetUtils.bitSetFrom(baseDeletionQualities[i]);
        }

        return new CovariateValues(mismatches, insertions, deletions);
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return Byte.parseByte(str);
    }

    @Override
    public String keyFromBitSet(BitSet key) {
        return String.format("%d", BitSetUtils.longFrom(key));
    }

    @Override
    public BitSet bitSetFromKey(Object key) {        
        return (key instanceof String) ? BitSetUtils.bitSetFrom(Byte.parseByte((String) key)) : BitSetUtils.bitSetFrom((Byte) key);
    }

    @Override
    public int numberOfBits() {
        return BitSetUtils.numberOfBitsToRepresent(QualityUtils.MAX_QUAL_SCORE);
    }
}
