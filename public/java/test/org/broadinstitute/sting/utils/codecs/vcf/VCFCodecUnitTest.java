/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// our package
package org.broadinstitute.sting.utils.codecs.vcf;


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class VCFCodecUnitTest extends BaseTest {

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class AlleleClippingTestProvider extends TestDataProvider {
        final String ref;
        final List<Allele> alleles = new ArrayList<Allele>();
        final int expectedClip;

        private AlleleClippingTestProvider(final int expectedClip, final String ref, final String ... alleles) {
            super(AlleleClippingTestProvider.class);
            this.ref = ref;
            for ( final String allele : alleles )
                this.alleles.add(Allele.create(allele));
            this.expectedClip = expectedClip;
        }

        @Override
        public String toString() {
            return String.format("ref=%s allele=%s reverse clip %d", ref, alleles, expectedClip);
        }
    }

    @DataProvider(name = "AlleleClippingTestProvider")
    public Object[][] MakeAlleleClippingTest() {
        // pair clipping
        new AlleleClippingTestProvider(0, "ATT", "CCG");
        new AlleleClippingTestProvider(1, "ATT", "CCT");
        new AlleleClippingTestProvider(2, "ATT", "CTT");
        new AlleleClippingTestProvider(2, "ATT", "ATT");  // cannot completely clip allele

        // triplets
        new AlleleClippingTestProvider(0, "ATT", "CTT", "CGG");
        new AlleleClippingTestProvider(1, "ATT", "CTT", "CGT"); // the T can go
        new AlleleClippingTestProvider(2, "ATT", "CTT", "CTT"); // both Ts can go

        return AlleleClippingTestProvider.getTests(AlleleClippingTestProvider.class);
    }


    @Test(dataProvider = "AlleleClippingTestProvider")
    public void TestAlleleClipping(AlleleClippingTestProvider cfg) {
        int result = AbstractVCFCodec.computeReverseClipping(cfg.alleles, cfg.ref, 0, 1);
        Assert.assertEquals(result, cfg.expectedClip);
    }
}