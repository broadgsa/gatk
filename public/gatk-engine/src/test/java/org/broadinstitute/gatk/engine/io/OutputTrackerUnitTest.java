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

package org.broadinstitute.gatk.engine.io;

import org.broadinstitute.gatk.engine.io.stubs.OutputStreamStub;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class OutputTrackerUnitTest extends BaseTest {

    private final OutputTracker tracker = new DirectOutputTracker();
    private File unwriteableDir = null;
    private File untraversableDir = null;

    @BeforeClass
    public void createDirectories() {
        unwriteableDir = new File("unwriteable");
        unwriteableDir.deleteOnExit();
        unwriteableDir.mkdir();
        unwriteableDir.setWritable(false);

        untraversableDir = new File("untraversable");
        untraversableDir.deleteOnExit();
        untraversableDir.mkdir();
        untraversableDir.setExecutable(false);
    }

    @DataProvider(name = "BadOutputPaths")
    public Object[][] makeBadOutputPaths() {
        return new Object[][] {new String[] {"thisDirectoryDoesNotExist/stub.txt"},
                new String[] {"/thisDirectoryDoesNotExist/dummy.txt"},
                new String[] {unwriteableDir.getAbsolutePath()+"/dummy.txt"},
                new String[] {untraversableDir.getAbsolutePath()+"/dummy.txt"}};
    }

    @DataProvider(name = "GoodOutputPaths")
    public Object[][] makeGoodOutputPaths() {
        return new Object[][] {new String[] {publicTestDir+"stub.txt"},
                new String[] {"dummy.txt"}};
    }

    @Test(dataProvider = "BadOutputPaths", expectedExceptions = UserException.CouldNotCreateOutputFile.class)
    public void testInvalidOutputPath(final String path) {
        tracker.validateOutputPath(new OutputStreamStub(new File(path)));
    }

    @Test(dataProvider = "GoodOutputPaths")
    public void testValidOutputPath(final String path) {
        tracker.validateOutputPath(new OutputStreamStub(new File(path)));
    }

    @Test
    public void testOutputPathWithNullFile() {
        tracker.validateOutputPath(new OutputStreamStub(System.out));
    }
}
