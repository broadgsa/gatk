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

package org.broadinstitute.gatk.tools.walkers.annotator;

/**
 * Created with IntelliJ IDEA.
 * User: farjoun
 * Date: 6/5/13
 * Time: 2:31 PM
 * To change this template use File | Settings | File Templates.
 */


import org.broadinstitute.gatk.tools.walkers.annotator.SnpEff.EffectType;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class SnpEffUtilUnitTest {


    @DataProvider(name="effects")
    public Object[][] childParentpairs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{EffectType.GENE,EffectType.CHROMOSOME});
        tests.add(new Object[]{EffectType.UTR_3_PRIME,EffectType.TRANSCRIPT});
        tests.add(new Object[]{EffectType.CODON_CHANGE,EffectType.CDS});
        tests.add(new Object[]{EffectType.STOP_GAINED,EffectType.EXON});
        tests.add(new Object[]{EffectType.SYNONYMOUS_START,EffectType.TRANSCRIPT});
        tests.add(new Object[]{EffectType.FRAME_SHIFT,EffectType.CDS});
        tests.add(new Object[]{EffectType.UPSTREAM,EffectType.INTERGENIC});
        tests.add(new Object[]{EffectType.SPLICE_SITE_DONOR,EffectType.INTRON});
        tests.add(new Object[]{EffectType.SPLICE_SITE_ACCEPTOR,EffectType.INTRON});
        tests.add(new Object[]{EffectType.STOP_LOST,EffectType.NON_SYNONYMOUS_CODING});
        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name="self")
    public Object[][] childEqualsParentpairs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for(EffectType type:EffectType.values()){
            tests.add(new Object[]{type,type});
        }
        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name="noneffects")
    public Object[][] nonchildParentpairs() {
        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{EffectType.START_GAINED,EffectType.NON_SYNONYMOUS_CODING});
        tests.add(new Object[]{EffectType.GENE,EffectType.NONE});
        tests.add(new Object[]{EffectType.UTR_3_PRIME,EffectType.CDS});
        tests.add(new Object[]{EffectType.CODON_CHANGE,EffectType.REGULATION});
        tests.add(new Object[]{EffectType.DOWNSTREAM,EffectType.REGULATION});
        tests.add(new Object[]{EffectType.SPLICE_SITE_ACCEPTOR,EffectType.EXON});
        tests.add(new Object[]{EffectType.START_GAINED,EffectType.SYNONYMOUS_START});
        tests.add(new Object[]{EffectType.NON_SYNONYMOUS_CODING,EffectType.DOWNSTREAM});
        tests.add(new Object[]{EffectType.CODON_DELETION,EffectType.INTRON});
        tests.add(new Object[]{EffectType.UTR_5_PRIME,EffectType.EXON_DELETED});
        tests.add(new Object[]{EffectType.INTRON,EffectType.NONE});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "effects")
    public void testSubType(EffectType subType,EffectType parentType) {
        Assert.assertTrue(SnpEffUtil.isSubTypeOf(subType,parentType),String.format("testing that %s is subtype of %s.",subType,parentType));
    }
    @Test(dataProvider = "self")
    public void testSubTypeSelf(EffectType subType,EffectType parentType) {
        Assert.assertTrue(SnpEffUtil.isSubTypeOf(subType,parentType),String.format("testing that %s is subtype of %s.",subType,parentType));
    }
    @Test(dataProvider = "effects")
    public void testNonSubTypeSelf(EffectType parentType,EffectType subType) {
        Assert.assertTrue(!SnpEffUtil.isSubTypeOf(subType,parentType),String.format("testing that %s is subtype of %s.",subType,parentType));
    }
    @Test(dataProvider = "noneffects")
    public void testNonSubType(EffectType subType,EffectType parentType) {
        Assert.assertTrue(!SnpEffUtil.isSubTypeOf(subType, parentType), String.format("testing that %s is NOT subtype of %s.", subType, parentType));
        Assert.assertTrue(!SnpEffUtil.isSubTypeOf(parentType,subType), String.format("testing that %s is NOT subtype of %s.", parentType,subType));
    }
}
