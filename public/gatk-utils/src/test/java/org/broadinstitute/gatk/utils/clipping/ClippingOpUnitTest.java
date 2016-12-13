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

package org.broadinstitute.gatk.utils.clipping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class ClippingOpUnitTest extends BaseTest {

    @Test (dataProvider = "SoftClipedReadsNewStart")
    public void testGetNewAlignmentStartOffset(final String preClip, final String postClip, final int expectedResult) {
        Cigar cpreClip = TextCigarCodec.decode(preClip);
        Cigar cpostClip = TextCigarCodec.decode(postClip);
        Assert.assertEquals(ClippingOp.getNewAlignmentStartOffset(cpostClip, cpreClip), expectedResult,
                "getNewAlignmentStartOffset returned "+ClippingOp.getNewAlignmentStartOffset(cpostClip, cpreClip)+
                        " when "+expectedResult+" was expected for "+preClip.toString()+" which was clipped to "+postClip.toString());
    }

    // provider fields: cigar string before clip, cigar string after clip, expected result for getNewAlignmentStartOffset() method
    @DataProvider(name = "SoftClipedReadsNewStart")
    public Object[][] makeRevertSoftClipsBeforeContig() {
        return new Object[][] {
                {"70M", "10S60M", 10},
                {"70M", "60M10S", 0},

                {"30M10N30M", "30S10N30M", 30},
                {"30M10N30M", "30S5N30M", 35},
                {"30M10N30M", "30M10N30S", 0},
                {"30M10N30M", "30S30M", 40},
                {"30M10N30M", "30M30S", 0},
                {"30M10N30M", "15S15M10N30M", 15},
                {"30M10N30M", "30M10N15M15S", 0},
                // Testing multiple sequential reference but not sequence consuming base types
                {"10N10D40M", "40M", 20},
                {"10N10D40M", "20S20M", 40},
                {"10N10D40M", "5N10D40M", 5},
                {"10N10D40M", "5D40M", 15},
                {"10N10D40M", "10N10D20M20S", 0},

                {"10S10N20M10N10M", "10S20M10N10M", 10},
                {"10S10N20M10N10M", "10S10N20M10S", 0},

                {"10S10I20M10I10M", "20S20M10I10M", 0},
                {"10S10I20M10I10M", "10S10I20M20S", 0},
                {"10S10I20M10I10M", "15S5I20M10I10M", 0},

                {"10H60M", "10H10S50M", 10},
                {"10H60M", "10H50M10S", 0},
                {"10H10S50M", "10H20S40M", 10},
                {"10X60M", "20S50M", 20},
                {"10I40N20M","10S20M",40}
        };
    }
}
