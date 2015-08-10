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
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class ArgumentMatchSourceUnitTest extends BaseTest {
    @Test
    public void testCommandLine() {
        ArgumentMatchSource source = ArgumentMatchSource.COMMAND_LINE;
        Assert.assertEquals(source.getType(), ArgumentMatchSourceType.CommandLine);
        Assert.assertNull(source.getDescription());
    }

    @Test
    public void testFile() {
        File f = new File("test");
        ArgumentMatchSource source = new ArgumentMatchFileSource(f);
        Assert.assertEquals(source.getType(), ArgumentMatchSourceType.Provider);
        Assert.assertEquals(source.getDescription(), "file " + f.getAbsolutePath());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullFile() {
        new ArgumentMatchSource(null);
    }

    @Test
    public void testEquals() {
        ArgumentMatchSource cmdLine = ArgumentMatchSource.COMMAND_LINE;
        ArgumentMatchSource fileA = new ArgumentMatchFileSource(new File("a"));
        ArgumentMatchSource fileB = new ArgumentMatchFileSource(new File("b"));

        Assert.assertFalse(cmdLine.equals(null));

        Assert.assertTrue(cmdLine.equals(cmdLine));
        Assert.assertFalse(cmdLine.equals(fileA));
        Assert.assertFalse(cmdLine.equals(fileB));

        Assert.assertFalse(fileA.equals(cmdLine));
        Assert.assertTrue(fileA.equals(fileA));
        Assert.assertFalse(fileA.equals(fileB));

        Assert.assertFalse(fileB.equals(cmdLine));
        Assert.assertFalse(fileB.equals(fileA));
        Assert.assertTrue(fileB.equals(fileB));
    }

    @Test
    public void testCompareTo() {
        ArgumentMatchSource cmdLine = ArgumentMatchSource.COMMAND_LINE;
        ArgumentMatchSource fileA = new ArgumentMatchFileSource(new File("a"));
        ArgumentMatchSource fileB = new ArgumentMatchFileSource(new File("b"));

        Assert.assertTrue(cmdLine.compareTo(cmdLine) == 0);
        Assert.assertTrue(cmdLine.compareTo(fileA) < 0);
        Assert.assertTrue(cmdLine.compareTo(fileB) < 0);

        Assert.assertTrue(fileA.compareTo(cmdLine) > 0);
        Assert.assertTrue(fileA.compareTo(fileA) == 0);
        Assert.assertTrue(fileA.compareTo(fileB) < 0);

        Assert.assertTrue(fileB.compareTo(cmdLine) > 0);
        Assert.assertTrue(fileB.compareTo(fileA) > 0);
        Assert.assertTrue(fileB.compareTo(fileB) == 0);
    }

    @Test(expectedExceptions = NullPointerException.class)
    public void testCompareToNull() {
        ArgumentMatchSource.COMMAND_LINE.compareTo(null);
    }
}
