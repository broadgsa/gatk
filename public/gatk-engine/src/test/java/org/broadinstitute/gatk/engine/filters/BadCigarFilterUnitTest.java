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

package org.broadinstitute.gatk.engine.filters;

import htsjdk.samtools.Cigar;
import org.broadinstitute.gatk.utils.clipping.ReadClipperTestUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.List;

/**
 * Checks that the Bad Cigar filter works for all kinds of wonky cigars
 *
 * @author Mauricio Carneiro
 * @since 3/20/12
 */
public class BadCigarFilterUnitTest {

    public static final String[] BAD_CIGAR_LIST = {
            "2D4M",               // starting with multiple deletions
            "4M2D",               // ending with multiple deletions
            "3M1I1D",             // adjacent indels AND ends in deletion
            "1M1I1D2M",           // adjacent indels I->D
            "1M1D2I1M",           // adjacent indels D->I
            "1M1I2M1D",           // ends in single deletion with insertion in the middle
            "4M1D",               // ends in single deletion
            "1D4M",               // starts with single deletion
            "2M1D1D2M",           // adjacent D's
            "1M1I1I1M",           // adjacent I's
            "1H1D4M",             // starting with deletion after H
            "1S1D3M",             // starting with deletion after S
            "1H1S1D3M",           // starting with deletion after HS
            "4M1D1H",             // ending with deletion before H
            "3M1D1S",             // ending with deletion before S
            "3M1D1S1H",           // ending with deletion before HS
            "10M2H10M",           // H in the middle
            "10M2S10M",           // S in the middle
            "1H1S10M2S10M1S1H",    // deceiving S in the middle
            "1H1S10M2H10M1S1H"    // deceiving H in the middle
    };

    BadCigarFilter filter;

    @BeforeClass
    public void init() {
        filter = new BadCigarFilter();
    }

    @Test(enabled = true)
    public void testWonkyCigars() {
        for (String cigarString : BAD_CIGAR_LIST) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigarString, 0);
            Assert.assertTrue(filter.filterOut(read), read.getCigarString());
        }
    }

    @Test(enabled = true)
    public void testReadCigarLengthMismatch() {
        GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar("4M", 1);
        Assert.assertTrue(filter.filterOut(read), read.getCigarString());
    }

    @Test(enabled = true)
    public void testGoodCigars() {
        List<Cigar> cigarList = ReadClipperTestUtils.generateCigarList(10);
        for (Cigar cigar : cigarList) {
            GATKSAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar, 0);
            Assert.assertFalse(filter.filterOut(read), read.getCigarString());
        }
    }
}
