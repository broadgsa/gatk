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

package org.broadinstitute.gatk.utils.collections;


// the imports for unit testing.

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;


/**
 * Basic unit test for DefaultHashMap
 */
public class DefaultHashMapUnitTest extends BaseTest {
    DefaultHashMap<String, Double> empty, hasOne, hasTen;
    Double initialDefault = 10.0;

    @BeforeMethod
    public void before() {
        empty = new DefaultHashMap<String, Double>(initialDefault);

        hasOne = new DefaultHashMap<String, Double>(initialDefault);
        hasOne.put("1", .1);

        hasTen = new DefaultHashMap<String, Double>(initialDefault);
        for (Integer i = 1; i <= 10; i++) {
            hasTen.put(i.toString(), i.doubleValue() / 10);
        }
    }

    @Test
    public void testBasicSizes() {
        logger.warn("Executing testBasicSizes");

        Assert.assertEquals(0, empty.size());
        Assert.assertEquals(1, hasOne.size());
        Assert.assertEquals(10, hasTen.size());
    }

    @Test
    public void testTenElements() {
        logger.warn("Executing testTenElements");

        for (Integer i = 1; i <= 10; i++) {
            Assert.assertEquals(i.doubleValue() / 10, hasTen.get(i.toString()));
        }
        Assert.assertEquals(initialDefault, hasTen.get("0"));
    }

    @Test
    public void testClear() {
        logger.warn("Executing testClear");

        empty.clear();
        hasOne.clear();
        hasTen.clear();

        Assert.assertEquals(0, empty.size());
        Assert.assertEquals(0, hasOne.size());
        Assert.assertEquals(0, hasTen.size());
    }


    @Test
    public void testSettingTenElements() {
        logger.warn("Executing testSettingTenElements");

        Assert.assertEquals(10, hasTen.size());
        for (Integer i = 1; i <= 10; i++) {
            hasTen.put(i.toString(), i.doubleValue());
        }

        Assert.assertEquals(10, hasTen.size());
        for (Integer i = 1; i <= 10; i++) {
            Assert.assertEquals(i.doubleValue(), hasTen.get(i.toString()));
        }
    }

    @Test
    public void testSettingDefault() {
        logger.warn("Executing testSettingDefault");

        Assert.assertEquals(initialDefault, empty.get("0"));
        Assert.assertEquals(initialDefault, hasOne.get("0"));
        Assert.assertEquals(initialDefault, hasTen.get("0"));

        empty.setDefaultValue(2 * initialDefault);
        hasOne.setDefaultValue(2 * initialDefault);
        hasTen.setDefaultValue(2 * initialDefault);

        Assert.assertEquals(2 * initialDefault, empty.get("0"));
        Assert.assertEquals(2 * initialDefault, hasOne.get("0"));
        Assert.assertEquals(2 * initialDefault, hasTen.get("0"));

    }

    @Test
    public void testAdd() {
        logger.warn("Executing testAdd");

        Assert.assertEquals(0, empty.size());

        Double x = 1.0;
        empty.put(x.toString(), x / 10);
        Assert.assertEquals(1, empty.size());
        Assert.assertEquals(.1, empty.get(x.toString()));

        x = 2.0;
        empty.put(x.toString(), x / 10);
        Assert.assertEquals(2, empty.size());
        Assert.assertEquals(.2, empty.get(x.toString()));

    }

    @Test
    public void testUnset() {
        logger.warn("Executing testUnset1");

        Assert.assertEquals(10, hasTen.size());
        Assert.assertEquals(.9, hasTen.get("9"));

        hasTen.remove("9");

        Assert.assertEquals(9, hasTen.size());
        Assert.assertEquals(initialDefault, hasTen.get("9"));

        hasTen.remove("1");

        Assert.assertEquals(8, hasTen.size());
        Assert.assertEquals(initialDefault, hasTen.get("1"));

    }
}
