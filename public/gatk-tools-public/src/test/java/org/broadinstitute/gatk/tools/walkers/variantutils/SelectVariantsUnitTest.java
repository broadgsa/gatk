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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class SelectVariantsUnitTest extends BaseTest {

    ///////////////////////////////////////////////////////////
    // Tests for maxIndelSize and minIndelSize functionality //
    ///////////////////////////////////////////////////////////

    @DataProvider(name = "MaxMinIndelSize")
    public Object[][] MaxMinIndelSizeTestData() {

        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int size : Arrays.asList(1, 3, 10, 100) ) {
            for ( final int otherSize : Arrays.asList(0, 1) ) {
                for ( final int max : Arrays.asList(0, 1, 5, 50, 100000) ) {
                    for ( final int min : Arrays.asList(0, 1, 5, 50) ) {
                        for (final String op : Arrays.asList("D", "I")) {
                            tests.add(new Object[]{size, otherSize, max, min, op});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MaxMinIndelSize")
    public void maxIndelSizeTest(final int size, final int otherSize, final int max, final int min, final String op) {

        final byte[] largerAllele = Utils.dupBytes((byte) 'A', size+1);
        final byte[] smallerAllele = Utils.dupBytes((byte) 'A', 1);

        final List<Allele> alleles = new ArrayList<Allele>(2);
        final Allele ref = Allele.create(op.equals("I") ? smallerAllele : largerAllele, true);
        final Allele alt = Allele.create(op.equals("D") ? smallerAllele : largerAllele, false);
        alleles.add(ref);
        alleles.add(alt);

        final VariantContext vc = new VariantContextBuilder("test", "1", 10, 10 + ref.length() - 1, alleles).make();

        boolean hasIndelTooLargeOrSmall = SelectVariants.containsIndelLargerOrSmallerThan(vc, max, min);
        Assert.assertEquals(hasIndelTooLargeOrSmall, size > max || size < min);
    }

}