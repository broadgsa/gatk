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


import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.EnumMap;


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

    @Test (expectedExceptions = UserException.MalformedVCF.class)
    public void testErrorBadFormat() {
        GenotypeLikelihoods gl = new GenotypeLikelihoods("adf,b,c");
        gl.getAsVector();
    }

    @Test
    public void testGetAsMap(){
        GenotypeLikelihoods gl = new GenotypeLikelihoods(v);
        //Log scale
        EnumMap<Genotype.Type,Double> glMap = gl.getAsMap(false);
        Assert.assertEquals(v[Genotype.Type.HOM_REF.ordinal()-1],glMap.get(Genotype.Type.HOM_REF));
        Assert.assertEquals(v[Genotype.Type.HET.ordinal()-1],glMap.get(Genotype.Type.HET));
        Assert.assertEquals(v[Genotype.Type.HOM_VAR.ordinal()-1],glMap.get(Genotype.Type.HOM_VAR));

        //Linear scale
        glMap = gl.getAsMap(true);
        double [] vl = MathUtils.normalizeFromLog10(v);
        Assert.assertEquals(vl[Genotype.Type.HOM_REF.ordinal()-1],glMap.get(Genotype.Type.HOM_REF));
        Assert.assertEquals(vl[Genotype.Type.HET.ordinal()-1],glMap.get(Genotype.Type.HET));
        Assert.assertEquals(vl[Genotype.Type.HOM_VAR.ordinal()-1],glMap.get(Genotype.Type.HOM_VAR));

        //Test missing likelihoods
        gl = new GenotypeLikelihoods(".");
        glMap = gl.getAsMap(false);
        Assert.assertNull(glMap);

    }

    @Test
    public void testCalculateNumLikelihoods() {    
        
        for (int nAlleles=2; nAlleles<=5; nAlleles++)
            // simplest case: diploid
            Assert.assertEquals(GenotypeLikelihoods.calculateNumLikelihoods(nAlleles, 2), nAlleles*(nAlleles+1)/2);

        // some special cases: ploidy = 20, #alleles = 4
        Assert.assertEquals(GenotypeLikelihoods.calculateNumLikelihoods(4, 20), 1771);
    }
    
    @Test
    public void testGetLog10GQ(){
        GenotypeLikelihoods gl = new GenotypeLikelihoods(vPLString);

        //GQ for the best guess genotype
        Assert.assertEquals(gl.getLog10GQ(Genotype.Type.HET),-3.9);

        double[] test = MathUtils.normalizeFromLog10(gl.getAsVector());

        //GQ for the other genotypes
        Assert.assertEquals(gl.getLog10GQ(Genotype.Type.HOM_REF), Math.log10(1.0 - test[Genotype.Type.HOM_REF.ordinal()-1]));
        Assert.assertEquals(gl.getLog10GQ(Genotype.Type.HOM_VAR), Math.log10(1.0 - test[Genotype.Type.HOM_VAR.ordinal()-1]));

       //Test missing likelihoods
        gl = new GenotypeLikelihoods(".");
        Assert.assertEquals(gl.getLog10GQ(Genotype.Type.HOM_REF),Double.NEGATIVE_INFINITY);
        Assert.assertEquals(gl.getLog10GQ(Genotype.Type.HET),Double.NEGATIVE_INFINITY);
        Assert.assertEquals(gl.getLog10GQ(Genotype.Type.HOM_VAR),Double.NEGATIVE_INFINITY);

    }

    @Test
    public void testgetQualFromLikelihoods() {
        double[] likelihoods = new double[]{-1, 0, -2};
        // qual values we expect for each possible "best" genotype
        double[] expectedQuals = new double[]{-0.04100161, -1, -0.003930294};

        for ( int i = 0; i < likelihoods.length; i++ ) {
            Assert.assertEquals(GenotypeLikelihoods.getQualFromLikelihoods(i, likelihoods), expectedQuals[i], 1e-6,
                    "GQ value for genotype " + i + " was not calculated correctly");
        }
    }

    private void assertDoubleArraysAreEqual(double[] v1, double[] v2) {
        Assert.assertEquals(v1.length, v2.length);
        for ( int i = 0; i < v1.length; i++ ) {
            Assert.assertEquals(v1[i], v2[i], 1e-6);
        }
    }

    @Test
    public void testCalculatePLindex(){
        int counter = 0;
        for ( int i = 0; i <= 3; i++ ) {
            for ( int j = i; j <= 3; j++ ) {
                Assert.assertEquals(GenotypeLikelihoods.calculatePLindex(i, j), GenotypeLikelihoods.PLindexConversion[counter++], "PL index of alleles " + i + "," + j + " was not calculated correctly");
            }
        }
    }

    @Test
    public void testGetAllelePair(){
        allelePairTest(0, 0, 0);
        allelePairTest(1, 0, 1);
        allelePairTest(2, 1, 1);
        allelePairTest(3, 0, 2);
        allelePairTest(4, 1, 2);
        allelePairTest(5, 2, 2);
        allelePairTest(6, 0, 3);
        allelePairTest(7, 1, 3);
        allelePairTest(8, 2, 3);
        allelePairTest(9, 3, 3);
    }
        
    private void allelePairTest(int PLindex, int allele1, int allele2) {
        Assert.assertEquals(GenotypeLikelihoods.getAllelePair(PLindex).alleleIndex1, allele1, "allele index " + allele1 + " from PL index " + PLindex + " was not calculated correctly");
        Assert.assertEquals(GenotypeLikelihoods.getAllelePair(PLindex).alleleIndex2, allele2, "allele index " + allele2 + " from PL index " + PLindex + " was not calculated correctly");
    }
}