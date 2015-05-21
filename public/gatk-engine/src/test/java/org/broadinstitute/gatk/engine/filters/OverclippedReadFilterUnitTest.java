/*
* Copyright (c) 2012 The Broad Institute
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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


/**
 * Tests for the OverclippedReadFilter
 */
public class OverclippedReadFilterUnitTest extends ReadFilterTest {

    @Test(enabled = true, dataProvider= "OverclippedDataProvider")
    public void testOverclippedFilter(final String cigarString, final boolean expectedResult) {

       final OverclippedReadFilter filter = new OverclippedReadFilter();
       final SAMRecord read = buildSAMRecord(cigarString);
       Assert.assertEquals(filter.filterOut(read), expectedResult, cigarString);
    }

    private SAMRecord buildSAMRecord(final String cigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        return this.createRead(cigar, 1, 0, 10);
    }

    @DataProvider(name= "OverclippedDataProvider")
    public Iterator<Object[]> overclippedDataProvider() {
        final List<Object[]> result = new LinkedList<Object[]>();

        result.add(new Object[] { "1S10M1S", true });
        result.add(new Object[] { "1S10X1S", true });
        result.add(new Object[] { "1H1S10M1S1H", true });
        result.add(new Object[] { "1S40M1S", false });
        result.add(new Object[] { "1S40X1S", false });
        result.add(new Object[] { "1H10M1S", false });
        result.add(new Object[] { "1S10M1H", false });
        result.add(new Object[] { "10M1S", false });
        result.add(new Object[] { "1S10M", false });
        result.add(new Object[] { "1S10M10D10M1S", true });
        result.add(new Object[] { "1S1M40I1S", false });
        result.add(new Object[] { "1S10I1S", true });
        result.add(new Object[] { "1S40I1S", false });

        return result.iterator();
    }
}
