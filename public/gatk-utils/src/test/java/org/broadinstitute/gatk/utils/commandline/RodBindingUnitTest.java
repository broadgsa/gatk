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

package org.broadinstitute.gatk.utils.commandline;

import org.broadinstitute.gatk.utils.BaseTest;
import htsjdk.variant.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.BeforeMethod;

/**
 * Test suite for the parsing engine.
 */
public class RodBindingUnitTest extends BaseTest {
    Tags mytags = new Tags();

    @BeforeMethod
    public void setUp() {
        RodBinding.resetNameCounter();
    }

    @Test
    public void testStandardRodBinding() {
        RodBinding<VariantContext> b = new RodBinding<VariantContext>(VariantContext.class, "b", "foo", "vcf", mytags);
        Assert.assertEquals(b.getName(), "b");
        Assert.assertEquals(b.getType(), VariantContext.class);
        Assert.assertEquals(b.getSource(), "foo");
        Assert.assertEquals(b.getTribbleType(), "vcf");
        Assert.assertEquals(b.isBound(), true);
    }

    @Test
    public void testUnboundRodBinding() {
        RodBinding<VariantContext> u = RodBinding.makeUnbound(VariantContext.class);
        Assert.assertEquals(u.getName(), RodBinding.UNBOUND_VARIABLE_NAME);
        Assert.assertEquals(u.getSource(), RodBinding.UNBOUND_SOURCE);
        Assert.assertEquals(u.getType(), VariantContext.class);
        Assert.assertEquals(u.getTribbleType(), RodBinding.UNBOUND_TRIBBLE_TYPE);
        Assert.assertEquals(u.isBound(), false);
    }

    @Test
    public void testMultipleBindings() {
        String name = "binding";
        RodBinding<VariantContext> b1 = new RodBinding<VariantContext>(VariantContext.class, name, "foo", "vcf", mytags);
        Assert.assertEquals(b1.getName(), name);
        Assert.assertEquals(b1.getType(), VariantContext.class);
        Assert.assertEquals(b1.getSource(), "foo");
        Assert.assertEquals(b1.getTribbleType(), "vcf");
        Assert.assertEquals(b1.isBound(), true);

        RodBinding<VariantContext> b2 = new RodBinding<VariantContext>(VariantContext.class, name, "foo", "vcf", mytags);
        Assert.assertEquals(b2.getName(), name + "2");
        Assert.assertEquals(b2.getType(), VariantContext.class);
        Assert.assertEquals(b2.getSource(), "foo");
        Assert.assertEquals(b2.getTribbleType(), "vcf");
        Assert.assertEquals(b2.isBound(), true);
    }
}
