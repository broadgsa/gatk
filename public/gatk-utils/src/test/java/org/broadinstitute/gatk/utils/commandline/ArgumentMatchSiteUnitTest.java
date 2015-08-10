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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class ArgumentMatchSiteUnitTest {
    @Test
    public void testCommandLine() {
        ArgumentMatchSite site = new ArgumentMatchSite(ArgumentMatchSource.COMMAND_LINE, 1);
        Assert.assertEquals(site.getSource(), ArgumentMatchSource.COMMAND_LINE);
        Assert.assertEquals(site.getIndex(), 1);
    }

    @Test
    public void testFile() {
        ArgumentMatchSource source = new ArgumentMatchFileSource(new File("test"));
        ArgumentMatchSite site = new ArgumentMatchSite(source, 1);
        Assert.assertEquals(site.getSource(), source);
        Assert.assertEquals(site.getIndex(), 1);
    }

    @Test
    public void testEquals() {
        ArgumentMatchSource cmdLine = ArgumentMatchSource.COMMAND_LINE;
        ArgumentMatchSite site1 = new ArgumentMatchSite(cmdLine, 1);
        ArgumentMatchSite site2 = new ArgumentMatchSite(cmdLine, 2);

        Assert.assertFalse(site1.equals(null));

        Assert.assertTrue(site1.equals(site1));
        Assert.assertFalse(site1.equals(site2));

        Assert.assertFalse(site2.equals(site1));
        Assert.assertTrue(site2.equals(site2));
    }

    @Test
    public void testCompareTo() {
        ArgumentMatchSource cmdLine = ArgumentMatchSource.COMMAND_LINE;
        ArgumentMatchSite site1 = new ArgumentMatchSite(cmdLine, 1);
        ArgumentMatchSite site2 = new ArgumentMatchSite(cmdLine, 2);

        Assert.assertTrue(site1.compareTo(site1) == 0);
        Assert.assertTrue(site1.compareTo(site2) < 0);
        Assert.assertTrue(site2.compareTo(site1) > 0);
        Assert.assertTrue(site2.compareTo(site2) == 0);
    }

    @Test(expectedExceptions = NullPointerException.class)
    public void testCompareToNull() {
        new ArgumentMatchSite(ArgumentMatchSource.COMMAND_LINE, 0).compareTo(null);
    }
}
