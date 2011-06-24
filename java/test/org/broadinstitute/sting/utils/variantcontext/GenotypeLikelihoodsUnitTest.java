/*
 * Copyright (c) 2010, The Broad Institute
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
package org.broadinstitute.sting.utils.variantcontext;


// the imports for unit testing.


import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;


/**
 * Basic unit test for Genotype likelihoods objects
 */
public class GenotypeLikelihoodsUnitTest {
    double [] v = new double[]{-10.5, -1.25, -5.11};
    final static String vGLString = "-10.50,-1.25,-5.11";
    final static String vPLString = "93,0,39";

    @Test
    public void testFromVector2() {
        GenotypeLikelihoods gl = new GenotypeLikelihoods(v);
        assertDoubleArraysAreEqual(gl.getAsVector(), v);
        Assert.assertEquals(gl.getAsString(), vPLString);
    }

    @Test
    public void testFromString1() {
        GenotypeLikelihoods gl = new GenotypeLikelihoods(vPLString);
        assertDoubleArraysAreEqual(gl.getAsVector(), new double[]{-9.3, 0, -3.9});
        Assert.assertEquals(gl.getAsString(), vPLString);
    }

    @Test
    public void testFromString2() {
        GenotypeLikelihoods gl = GenotypeLikelihoods.fromGLField(vGLString);
        assertDoubleArraysAreEqual(gl.getAsVector(), v);
        Assert.assertEquals(gl.getAsString(), vPLString);
    }

    @Test (expectedExceptions = NumberFormatException.class)
    public void testErrorBadFormat() {
        GenotypeLikelihoods gl = new GenotypeLikelihoods("adf,b,c");
        gl.getAsVector();
    }

    private void assertDoubleArraysAreEqual(double[] v1, double[] v2) {
        Assert.assertEquals(v1.length, v2.length);
        for ( int i = 0; i < v1.length; i++ ) {
            Assert.assertEquals(v1[i], v2[i], 1e-6);
        }
    }
}