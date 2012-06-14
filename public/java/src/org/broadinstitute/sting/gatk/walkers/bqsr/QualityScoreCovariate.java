package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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

        long[] mismatches = new long[readLength];
        long[] insertions = new long[readLength];
        long[] deletions = new long[readLength];

        byte[] baseQualities = read.getBaseQualities();
        byte[] baseInsertionQualities = read.getBaseInsertionQualities();
        byte[] baseDeletionQualities = read.getBaseDeletionQualities();

        for (int i = 0; i < baseQualities.length; i++) {
            mismatches[i] = (long)baseQualities[i];
            insertions[i] = (long)baseInsertionQualities[i];
            deletions[i] = (long)baseDeletionQualities[i];
        }

        return new CovariateValues(mismatches, insertions, deletions);
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return Byte.parseByte(str);
    }

    @Override
    public String formatKey(final long key) {
        return String.format("%d", key);
    }

    @Override
    public long longFromKey(final Object key) {
        return (key instanceof String) ? (long)Byte.parseByte((String) key) : (long)(Byte) key;
    }

    @Override
    public int numberOfBits() {
        return BQSRKeyManager.numberOfBitsToRepresent(QualityUtils.MAX_QUAL_SCORE);
    }
}
