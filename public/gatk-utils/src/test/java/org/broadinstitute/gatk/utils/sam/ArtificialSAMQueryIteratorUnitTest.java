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

package org.broadinstitute.gatk.utils.sam;

import org.broadinstitute.gatk.utils.BaseTest;
import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;
import htsjdk.samtools.SAMRecord;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class ArtificialSAMQueryIteratorUnitTest
 *         <p/>
 *         a test for the ArtificialSAMQueryIterator class.
 */
public class ArtificialSAMQueryIteratorUnitTest extends BaseTest {

    @Test
    public void testWholeChromosomeQuery() {
        ArtificialSAMQueryIterator iter = ArtificialSAMUtils.queryReadIterator(1, 2, 100);
        iter.queryContained("chr1", 1, -1);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            count++;
        }
        assertEquals(count, 100);

    }

    @Test
    public void testContainedQueryStart() {
        ArtificialSAMQueryIterator iter = ArtificialSAMUtils.queryReadIterator(1, 2, 100);
        iter.queryContained("chr1", 1, 50);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            count++;
        }
        assertEquals(count, 1);

    }

    @Test
    public void testOverlappingQueryStart() {
        ArtificialSAMQueryIterator iter = ArtificialSAMUtils.queryReadIterator(1, 2, 100);
        iter.queryOverlapping("chr1", 1, 50);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            count++;
        }
        assertEquals(count, 50);

    }

    @Test
    public void testContainedQueryMiddle() {
        ArtificialSAMQueryIterator iter = ArtificialSAMUtils.queryReadIterator(1, 2, 100);
        iter.queryContained("chr1", 25, 74);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            count++;
        }
        assertEquals(count, 1);

    }

    @Test
    public void testOverlappingQueryMiddle() {
        ArtificialSAMQueryIterator iter = ArtificialSAMUtils.queryReadIterator(1, 2, 100);
        iter.queryOverlapping("chr1", 25, 74);
        int count = 0;
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            count++;
        }
        assertEquals(count, 50);

    }

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testUnknownChromosome() {
        ArtificialSAMQueryIterator iter = ArtificialSAMUtils.queryReadIterator(1, 2, 100);
        iter.queryOverlapping("chr621", 25, 74);         
    }
}
