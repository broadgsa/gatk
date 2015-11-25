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

package org.broadinstitute.gatk.engine.datasources.reads;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class SAMReaderIDUnitTest extends BaseTest {

    @Test
    public void testSAMReaderIDHashingAndEquality() {
        // Test to make sure that two SAMReaderIDs that point at the same file via an absolute vs. relative
        // path are equal according to equals() and have the same hash code
        final File relativePathToBAMFile = new File(publicTestDir + "exampleBAM.bam");
        final File absolutePathToBAMFile = new File(relativePathToBAMFile.getAbsolutePath());
        final SAMReaderID relativePathSAMReaderID = new SAMReaderID(relativePathToBAMFile, new Tags());
        final SAMReaderID absolutePathSAMReaderID = new SAMReaderID(absolutePathToBAMFile, new Tags());

        Assert.assertEquals(relativePathSAMReaderID, absolutePathSAMReaderID, "Absolute-path and relative-path SAMReaderIDs not equal according to equals()");
        Assert.assertEquals(relativePathSAMReaderID.hashCode(), absolutePathSAMReaderID.hashCode(), "Absolute-path and relative-path SAMReaderIDs have different hash codes");
    }
}
