/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import org.junit.Test;
import org.broadinstitute.sting.BaseTest;
import junit.framework.Assert;

/**
 * Testing framework for general purpose utilities class.
 *
 * @author hanna
 * @version 0.1
 */

public class UtilsTest extends BaseTest {

    @Test
    public void testDupStringNoChars() {
        String duped = Utils.dupString('a',0);
        Assert.assertEquals("dupString did not produce zero-length string", 0, duped.length());
    }

    @Test
    public void testDupStringOneChar() {
        String duped = Utils.dupString('b',1);
        Assert.assertEquals("dupString did not produce single character string", 1, duped.length());
        Assert.assertEquals("dupString character was incorrect", 'b', duped.charAt(0));
    }

    @Test
    public void testDupStringMultiChar() {
        String duped = Utils.dupString('c',5);
        Assert.assertEquals("dupString did not produce five character string", 5, duped.length());
        Assert.assertEquals("dupString string was incorrect", "ccccc", duped);
    }

}
