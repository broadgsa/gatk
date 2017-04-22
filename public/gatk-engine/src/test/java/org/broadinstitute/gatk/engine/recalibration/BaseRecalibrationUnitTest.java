/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.recalibration;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BaseRecalibrationUnitTest {

    @Test
    public void repeatedAndUnorderedFixedQualities() {
        // Test both repeated quals, and quals that aren't input in order
        List<Integer> quantizedQualsOrdered = Arrays.asList(11, 19);
        List<Integer> quantizedQualsUnordered = Arrays.asList(19, 11, 19, 19);

        // Unordered and Ordered qmapping should be identical
        byte[] qmappingUnordered = BaseRecalibration.constructStaticQuantizedMapping(quantizedQualsUnordered, true);
        byte[] qmappingOrdered = BaseRecalibration.constructStaticQuantizedMapping(quantizedQualsOrdered, true);
        Assert.assertEquals(qmappingOrdered.length, qmappingUnordered.length);
        for(int i = 0 ; i < qmappingUnordered.length ; i++) {
            Assert.assertEquals(qmappingOrdered[i], qmappingUnordered[i]);
        }
    }

    @Test
    public void nearestVsRoundDown() {
        List<Integer> fixedQuantizedQuals = Arrays.asList(10, 20, 30);

        byte[] qmappingRoundDown = BaseRecalibration.constructStaticQuantizedMapping(fixedQuantizedQuals, true);
        byte[] qmappingRoundNearest = BaseRecalibration.constructStaticQuantizedMapping(fixedQuantizedQuals, false);

        // Depending on rounding strategy, bin 19 should round to 10 or 20
        Assert.assertEquals(qmappingRoundDown[19], 10);
        Assert.assertEquals(qmappingRoundNearest[19], 20);

        // Regarless of rounding strategy, bin 21 should always round down to 20
        Assert.assertEquals(qmappingRoundDown[21], 20);
        Assert.assertEquals(qmappingRoundNearest[21], 20);
    }

    @Test
    public void onlyOneFixedQualUsed() {
        // Set all qualities to singleQual value (except for those below MIN_USABLE_Q_SCORE)
        int singleQual = 10;
        List<Integer> fixedQuantizedQuals = Arrays.asList(singleQual);

        byte[] qmapping = BaseRecalibration.constructStaticQuantizedMapping(fixedQuantizedQuals, true);

        for(int i = 0 ; i < qmapping.length ; i++) {
            if(i >= QualityUtils.MIN_USABLE_Q_SCORE) {
                Assert.assertEquals(qmapping[i], singleQual);
            }
            else {
                // Make sure that all values less than MIN_USABLE_Q_SCORE are preserved
                Assert.assertEquals(qmapping[i], i);
            }
        }
    }
}
