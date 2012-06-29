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
package org.broadinstitute.sting.utils.codecs.vcf;

import com.google.java.contract.Requires;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class VCFUtilsUnitTest extends BaseTest {
    // --------------------------------------------------------------------------------
    //
    // Test allele clipping
    //
    // --------------------------------------------------------------------------------

    private class ClipAllelesTest extends TestDataProvider {
        final int position;
        final int stop;
        final String ref;
        List<Allele> inputs;
        List<Allele> expected;

        @Requires("arg.length % 2 == 0")
        private ClipAllelesTest(final int position, final int stop, final String ... arg) {
            super(ClipAllelesTest.class);
            this.position = position;
            this.stop = stop;
            this.ref = arg[0];

            int n = arg.length / 2;
            inputs = new ArrayList<Allele>(n);
            expected = new ArrayList<Allele>(n);

            for ( int i = 0; i < n; i++ ) {
                final boolean ref = i % n == 0;
                inputs.add(Allele.create(arg[i], ref));
            }
            for ( int i = n; i < arg.length; i++ ) {
                final boolean ref = i % n == 0;
                expected.add(Allele.create(arg[i], ref));
            }
        }

        public String toString() {
            return String.format("ClipAllelesTest input=%s expected=%s", inputs, expected);
        }
    }
    @DataProvider(name = "ClipAllelesTest")
    public Object[][] makeClipAllelesTest() {
        // do no harm
        new ClipAllelesTest(10, 10, "A", "A");
        new ClipAllelesTest(10, 10, "A", "C", "A", "C");
        new ClipAllelesTest(10, 10, "A", "C", "G", "A", "C", "G");

        // insertions
        new ClipAllelesTest(10, 10, "A", "AA", "-", "A");
        new ClipAllelesTest(10, 10, "A", "AAA", "-", "AA");
        new ClipAllelesTest(10, 10, "A", "AG", "-", "G");

        // deletions
        new ClipAllelesTest(10, 11, "AA",  "A", "A",  "-");
        new ClipAllelesTest(10, 12, "AAA", "A", "AA", "-");
        new ClipAllelesTest(10, 11, "AG",  "A", "G",  "-");
        new ClipAllelesTest(10, 12, "AGG", "A", "GG", "-");

        // multi-allelic insertion and deletions
        new ClipAllelesTest(10, 11, "AA",  "A", "AAA", "A",  "-", "AA");
        new ClipAllelesTest(10, 11, "AA",  "A", "AAG", "A",  "-", "AG");
        new ClipAllelesTest(10, 10, "A",  "AA", "AAA", "-",  "A", "AA");
        new ClipAllelesTest(10, 10, "A",  "AA", "ACA", "-",  "A", "CA");
        new ClipAllelesTest(10, 12, "ACG", "ATC", "AGG", "CG",  "TC", "GG");
        new ClipAllelesTest(10, 11, "AC", "AT", "AG", "C",  "T", "G");

        // cannot be clipped
        new ClipAllelesTest(10, 11, "AC", "CT", "AG", "AC",  "CT", "AG");
        new ClipAllelesTest(10, 11, "AC", "CT", "GG", "AC",  "CT", "GG");

        // symbolic
        new ClipAllelesTest(10, 10, "A", "<DEL>", "A", "<DEL>");
        new ClipAllelesTest(50, 50, "G", "G]22:60]", "A", "G]22:60]");
        new ClipAllelesTest(51, 51, "T", "]22:55]T", "A", "]22:55]T");
        new ClipAllelesTest(52, 52, "C", "C[22:51[", "A", "C[22:51[");
        new ClipAllelesTest(60, 60, "A", "A]22:50]", "A", "A]22:50]");

        // complex substitutions
        new ClipAllelesTest(10, 10, "A", "GA", "A", "GA");

        return ClipAllelesTest.getTests(ClipAllelesTest.class);
    }

    @Test(dataProvider = "ClipAllelesTest")
    public void testClipAllelesTest(ClipAllelesTest cfg) {
        final VCFAlleleClipper.ClippedAlleles clipped = VCFAlleleClipper.clipAlleles(cfg.position, cfg.ref, cfg.inputs);
        Assert.assertNull(clipped.getError(), "Unexpected error occurred");
        Assert.assertEquals(clipped.getStop(), cfg.stop, "Clipped alleles stop");
        Assert.assertEquals(clipped.getClippedAlleles(), cfg.expected, "Clipped alleles");
    }


    // --------------------------------------------------------------------------------
    //
    // basic allele clipping test
    //
    // --------------------------------------------------------------------------------

    private class ReverseClippingPositionTestProvider extends TestDataProvider {
        final String ref;
        final List<Allele> alleles = new ArrayList<Allele>();
        final int expectedClip;

        private ReverseClippingPositionTestProvider(final int expectedClip, final String ref, final String... alleles) {
            super(ReverseClippingPositionTestProvider.class);
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

    @DataProvider(name = "ReverseClippingPositionTestProvider")
    public Object[][] makeReverseClippingPositionTestProvider() {
        // pair clipping
        new ReverseClippingPositionTestProvider(0, "ATT", "CCG");
        new ReverseClippingPositionTestProvider(1, "ATT", "CCT");
        new ReverseClippingPositionTestProvider(2, "ATT", "CTT");
        new ReverseClippingPositionTestProvider(2, "ATT", "ATT");  // cannot completely clip allele

        // triplets
        new ReverseClippingPositionTestProvider(0, "ATT", "CTT", "CGG");
        new ReverseClippingPositionTestProvider(1, "ATT", "CTT", "CGT"); // the T can go
        new ReverseClippingPositionTestProvider(2, "ATT", "CTT", "CTT"); // both Ts can go

        return ReverseClippingPositionTestProvider.getTests(ReverseClippingPositionTestProvider.class);
    }


    @Test(dataProvider = "ReverseClippingPositionTestProvider")
    public void testReverseClippingPositionTestProvider(ReverseClippingPositionTestProvider cfg) {
        int result = VCFAlleleClipper.computeReverseClipping(cfg.alleles, cfg.ref.getBytes(), 0, false);
        Assert.assertEquals(result, cfg.expectedClip);
    }
}
