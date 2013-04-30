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

package org.broadinstitute.sting.gatk.filters;

import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.Test;


/**
 * Tests for the MalformedReadFilter
 *
 * @author Eric Banks
 * @since 3/14/13
 */
public class MalformedReadFilterUnitTest {

    //////////////////////////////////////
    // Test the checkSeqStored() method //
    //////////////////////////////////////

    @Test(enabled = true)
    public void testcheckSeqStored () {

        final GATKSAMRecord goodRead = ArtificialSAMUtils.createArtificialRead(new byte[]{(byte)'A'}, new byte[]{(byte)'A'}, "1M");
        final GATKSAMRecord badRead = ArtificialSAMUtils.createArtificialRead(new byte[]{}, new byte[]{}, "1M");
        badRead.setReadString("*");

        Assert.assertTrue(MalformedReadFilter.checkSeqStored(goodRead, true));
        Assert.assertFalse(MalformedReadFilter.checkSeqStored(badRead, true));

        try {
            MalformedReadFilter.checkSeqStored(badRead, false);
            Assert.assertTrue(false, "We should have exceptioned out in the previous line");
        } catch (UserException e) { }
    }
}
