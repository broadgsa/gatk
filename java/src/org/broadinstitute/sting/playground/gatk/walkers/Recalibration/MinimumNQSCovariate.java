package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.QualityUtils;

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
 * Date: Nov 4, 2009
 */

public class MinimumNQSCovariate implements Covariate {

    public final static String ORIGINAL_QUAL_ATTRIBUTE_TAG = "OQ";
    protected boolean USE_ORIGINAL_QUALS;

    public MinimumNQSCovariate() { // empty constructor is required to instantiate covariate in CovariateCounterWalker and TableRecalibrationWalker
        USE_ORIGINAL_QUALS = false;
    }

    public MinimumNQSCovariate(boolean originalQuals) {
        USE_ORIGINAL_QUALS = originalQuals;
    }

    public Comparable getValue(SAMRecord read, int offset, char[] refBases) {
        byte[] quals = read.getBaseQualities();
        if ( USE_ORIGINAL_QUALS && read.getAttribute(ORIGINAL_QUAL_ATTRIBUTE_TAG) != null ) {
            Object obj = read.getAttribute(ORIGINAL_QUAL_ATTRIBUTE_TAG);
            if ( obj instanceof String )
                quals = QualityUtils.fastqToPhred((String)obj);
            else {
                throw new RuntimeException(String.format("Value encoded by %s in %s isn't a string!", ORIGINAL_QUAL_ATTRIBUTE_TAG, read.getReadName()));
            }
        }

        Integer minQual = (int)(quals[0]);
        for ( int qual : quals ) {
            if( qual < minQual ) {
                minQual = qual;
            }
        }
        return minQual;
    }
    
    public Comparable getValue(String str) {
        return Integer.parseInt( str );
    }

    public int estimatedNumberOfBins() {
        return 40;
    }

    public String toString() {
        return "Minimum Neighborhood Quality Score";
    }
}
