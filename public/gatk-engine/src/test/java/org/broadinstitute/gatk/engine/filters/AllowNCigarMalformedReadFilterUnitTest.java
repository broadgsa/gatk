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


import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;


/**
 * Tests for the {@link MalformedReadFilter} when the unsafe flag
 * {@link ValidationExclusion.TYPE#ALLOW_N_CIGAR_READS} is set.
 *
 * @author Valentin Ruano-Rubio
 * @since 6/6/13
 */
public class AllowNCigarMalformedReadFilterUnitTest extends MalformedReadFilterUnitTest {


    @Override
    protected ValidationExclusion composeValidationExclusion() {
        return new ValidationExclusion(Collections.singletonList(ValidationExclusion.TYPE.ALLOW_N_CIGAR_READS));
    }


    @Test(enabled = true,
            dataProvider= "UnsupportedCigarOperatorDataProvider")
    @CigarOperatorTest(CigarOperatorTest.Outcome.IGNORE)
    public void testCigarNOperatorFilterIgnore(final String cigarString) {

        final MalformedReadFilter filter = buildMalformedReadFilter(false);
        final SAMRecord nContainingCigarRead = buildSAMRecord(cigarString);
        Assert.assertFalse(filter.filterOut(nContainingCigarRead),
                "filters out N containing Cigar when it should ignore the fact");
    }

    @Test(enabled = false)
    @Override
    public void testCigarNOperatorFilterException(final String cigarString) {
        // Nothing to do here.
        // Just deactivates the parents test case.
    }







}
