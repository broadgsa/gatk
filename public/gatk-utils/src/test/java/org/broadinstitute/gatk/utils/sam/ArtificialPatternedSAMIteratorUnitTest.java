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
import static org.testng.Assert.assertTrue;
import static org.testng.Assert.fail;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;
import htsjdk.samtools.SAMFileHeader;
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
 * 
 * @author aaron 
 * 
 * Class ArtificialPatternedSAMIteratorUnitTest
 *
 * tests ArtificialPatternedSAMIterator, making sure that if you specify in order
 * you get reads in order, and if you specify out of order you get them out of order.  
 */
public class ArtificialPatternedSAMIteratorUnitTest extends BaseTest {

    // our artifical patterned iterator
    ArtificialPatternedSAMIterator iter;

    private int startingChr = 1;
    private int endingChr = 2;
    private int readCount = 100;
    private int DEFAULT_READ_LENGTH = ArtificialSAMUtils.DEFAULT_READ_LENGTH;
    SAMFileHeader header;

    @BeforeMethod
    public void before() {
        header = ArtificialSAMUtils.createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

    }
    @Test
    public void testInOrder() {
        iter = new ArtificialPatternedSAMIterator(startingChr,endingChr,readCount,0,header, ArtificialPatternedSAMIterator.PATTERN.IN_ORDER_READS);
        if (!iter.hasNext()) {
            fail("no reads in the ArtificialPatternedSAMIterator");
        }
        SAMRecord last = iter.next();
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            if (!(rec.getReferenceIndex() > last.getReferenceIndex()) && (rec.getAlignmentStart() <= last.getAlignmentStart())) {
                fail("read " + rec.getReadName() + " out of order compared to last read, " + last.getReadName());
            }
            last = rec;
        }

    }
    @Test
    public void testOutOfOrder() {
        int outOfOrderCount = 0;
        iter = new ArtificialPatternedSAMIterator(startingChr,endingChr,readCount,0,header, ArtificialPatternedSAMIterator.PATTERN.RANDOM_READS);
        if (!iter.hasNext()) {
            fail("no reads in the ArtificialPatternedSAMIterator");
        }
        SAMRecord last = iter.next();
        while (iter.hasNext()) {
            SAMRecord rec = iter.next();
            if (!(rec.getReferenceIndex() > last.getReferenceIndex()) && (rec.getAlignmentStart() <= last.getAlignmentStart())) {
                ++outOfOrderCount;
            }
            last = rec;
        }
        assertTrue(outOfOrderCount > 0);
    }


}
