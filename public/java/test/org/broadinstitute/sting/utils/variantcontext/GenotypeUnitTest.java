/*
 * Copyright (c) 2011, The Broad Institute
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


import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;


public class GenotypeUnitTest extends BaseTest {
    Allele A, Aref, T;

    @BeforeSuite
    public void before() {
        A = Allele.create("A");
        Aref = Allele.create("A", true);
        T = Allele.create("T");
    }

//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased) {
//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased, double[] log10Likelihoods) {
//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, double[] log10Likelihoods)
//    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError)
//    public Genotype(String sampleName, List<Allele> alleles)
//    public List<Allele> getAlleles()
//    public List<Allele> getAlleles(Allele allele)
//    public Allele getAllele(int i)
//    public boolean isPhased()
//    public int getPloidy()
//    public Type getType()
//    public boolean isHom()
//    public boolean isHomRef()
//    public boolean isHomVar()
//    public boolean isHet()
//    public boolean isNoCall()
//    public boolean isCalled()
//    public boolean isAvailable()
//    public boolean hasLikelihoods()
//    public GenotypeLikelihoods getLikelihoods()
//    public boolean sameGenotype(Genotype other)
//    public boolean sameGenotype(Genotype other, boolean ignorePhase)
//    public String getSampleName()
//    public boolean hasNegLog10PError()
//    public double getNegLog10PError()
//    public double getPhredScaledQual()
//    public boolean hasAttribute(String key)
//    public Object getAttribute(String key)
//    public Object getAttribute(String key, Object defaultValue)
//    public String getAttributeAsString(String key, String defaultValue)
//    public int getAttributeAsInt(String key, int defaultValue)
//    public double getAttributeAsDouble(String key, double  defaultValue)
//    public boolean getAttributeAsBoolean(String key, boolean  defaultValue)
}
