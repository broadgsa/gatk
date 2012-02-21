/*
 * Copyright (c) 2011 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 9/26/11
 */

public class ContextCovariate implements StandardCovariate {

    private int mismatchesContextSize;
    private int insertionsContextSize;
    private int  deletionsContextSize;

    private String mismatchesNoContext = "";
    private String insertionsNoContext = "";
    private String  deletionsNoContext = "";
    
    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        insertionsContextSize = RAC.INSERTIONS_CONTEXT_SIZE;
        deletionsContextSize = RAC.DELETIONS_CONTEXT_SIZE;

        if (mismatchesContextSize <= 0 || insertionsContextSize <= 0 || deletionsContextSize <= 0)
            throw new UserException(String.format("Context Size must be positive, if you don't want to use the context covariate, just turn it off instead. Mismatches: %d Insertions: %d Deletions:%d", mismatchesContextSize, insertionsContextSize, deletionsContextSize));

        // initialize no context strings given the size of the context for each covariate type
        mismatchesNoContext = makeAllNStringWithLength(mismatchesContextSize);
        insertionsNoContext = makeAllNStringWithLength(insertionsContextSize);
        deletionsNoContext  = makeAllNStringWithLength( deletionsContextSize);        
    }

    @Override
    public CovariateValues getValues(final GATKSAMRecord read) {
        int l = read.getReadLength();
        String[] mismatches = new String [l];
        String[] insertions = new String [l];
        String[]  deletions = new String [l];

        final boolean negativeStrand = read.getReadNegativeStrandFlag();
        byte[] bases = read.getReadBases();
        if (negativeStrand) {
            bases = BaseUtils.simpleReverseComplement(bases); //this is NOT in-place
        }
        for (int i = 0; i < read.getReadLength(); i++) {
            mismatches[i] = contextWith(bases, i, mismatchesContextSize, mismatchesNoContext);
            insertions[i] = contextWith(bases, i, insertionsContextSize, insertionsNoContext);
            deletions[i]  = contextWith(bases, i,  deletionsContextSize,  deletionsNoContext);
        }
        if (negativeStrand) {
            reverse(mismatches);
            reverse(insertions);
            reverse(deletions);
        }
        return new CovariateValues(mismatches, insertions, deletions);
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Comparable getValue(final String str) {
        return str;
    }

    /**
     * calculates the context of a base independent of the covariate mode
     *
     * @param bases           the bases in the read to build the context from
     * @param offset          the position in the read to calculate the context for
     * @param contextSize     context size to use building the context
     * @param noContextString string to return if the position is not far enough in the read to have a full context before.
     * @return
     */
    private String contextWith(byte [] bases, int offset, int contextSize, String noContextString) {
        return (offset < contextSize) ? noContextString : new String(Arrays.copyOfRange(bases, offset - contextSize, offset));
    } 
    
    private String makeAllNStringWithLength(int length) {
        String s = "";
        for (int i=0; i<length; i++)
            s += "N";
        return s;
    }

    /**
     * Reverses the given array in place.
     *
     * @param array any array
     */
    private static void reverse(final Comparable[] array) {
        final int arrayLength = array.length;
        for (int l = 0, r = arrayLength - 1; l < r; l++, r--) {
            final Comparable temp = array[l];
            array[l] = array[r];
            array[r] = temp;
        }
    }
}
