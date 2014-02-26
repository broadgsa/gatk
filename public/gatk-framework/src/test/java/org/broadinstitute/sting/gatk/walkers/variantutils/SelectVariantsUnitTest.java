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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class SelectVariantsUnitTest extends BaseTest {

    //////////////////////////////////////////
    // Tests for maxIndelSize functionality //
    //////////////////////////////////////////

    @DataProvider(name = "MaxIndelSize")
    public Object[][] MaxIndelSizeTestData() {

        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int size : Arrays.asList(1, 3, 10, 100) ) {
            for ( final int otherSize : Arrays.asList(0, 1) ) {
                for ( final int max : Arrays.asList(0, 1, 5, 50, 100000) ) {
                    for ( final String op : Arrays.asList("D", "I") ) {
                        tests.add(new Object[]{size, otherSize, max, op});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MaxIndelSize")
    public void maxIndelSizeTest(final int size, final int otherSize, final int max, final String op) {

        final byte[] largerAllele = Utils.dupBytes((byte) 'A', size+1);
        final byte[] smallerAllele = Utils.dupBytes((byte) 'A', 1);

        final List<Allele> alleles = new ArrayList<Allele>(2);
        final Allele ref = Allele.create(op.equals("I") ? smallerAllele : largerAllele, true);
        final Allele alt = Allele.create(op.equals("D") ? smallerAllele : largerAllele, false);
        alleles.add(ref);
        alleles.add(alt);
        if ( otherSize > 0 && otherSize != size ) {
            final Allele otherAlt = Allele.create(op.equals("D") ? Utils.dupBytes((byte) 'A', size-otherSize+1) : Utils.dupBytes((byte) 'A', otherSize+1), false);
            alleles.add(otherAlt);
        }

        final VariantContext vc = new VariantContextBuilder("test", "1", 10, 10 + ref.length() - 1, alleles).make();

        boolean hasTooLargeIndel = SelectVariants.containsIndelLargerThan(vc, max);
        Assert.assertEquals(hasTooLargeIndel, size > max);
    }

}