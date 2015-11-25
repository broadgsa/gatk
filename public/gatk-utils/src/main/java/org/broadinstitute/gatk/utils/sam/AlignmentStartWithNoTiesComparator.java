/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.sam;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import htsjdk.samtools.SAMRecord;

import java.util.Comparator;

public class AlignmentStartWithNoTiesComparator implements Comparator<SAMRecord> {
    @Requires("c1 >= 0 && c2 >= 0")
    @Ensures("result == 0 || result == 1 || result == -1")
    private int compareContigs(int c1, int c2) {
        if (c1 == c2)
            return 0;
        else if (c1 > c2)
            return 1;
        return -1;
    }

    @Requires("r1 != null && r2 != null")
    @Ensures("result == 0 || result == 1 || result == -1")
    public int compare(SAMRecord r1, SAMRecord r2) {
        int result;

        if (r1 == r2)
            result = 0;

        else if (r1.getReadUnmappedFlag())
            result = 1;
        else if (r2.getReadUnmappedFlag())
            result = -1;
        else {
            final int cmpContig = compareContigs(r1.getReferenceIndex(), r2.getReferenceIndex());

            if (cmpContig != 0)
                result = cmpContig;

            else {
                if (r1.getAlignmentStart() < r2.getAlignmentStart())
                    result = -1;
                else
                    result = 1;
            }
        }

        return result;
    }
}