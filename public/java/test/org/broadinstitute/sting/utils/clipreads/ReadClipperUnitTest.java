/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils.clipreads;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 9/28/11
 * Time: 9:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadClipperUnitTest extends BaseTest {

    // TODO: Add error messages on failed tests

    GATKSAMRecord read, expected;
    ReadClipper readClipper;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?"; //ASCII values = 33,43,53,63

    @BeforeClass
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, BASES.length());
        read.setReadUnmappedFlag(true);
        read.setReadBases(new String(BASES).getBytes());
        read.setBaseQualityString(new String(QUALS));

        readClipper = new ReadClipper(read);
    }

    @Test ( enabled = false )
    public void testHardClipBothEndsByReferenceCoordinates() {
        logger.warn("Executing testHardClipBothEndsByReferenceCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipBothEndsByReferenceCoordinates(0,0), new GATKSAMRecord(read.getHeader()));
        //clip 1 base
        expected = readClipper.hardClipBothEndsByReferenceCoordinates(0,3);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,3).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,3));
        Assert.assertEquals(expected.getCigarString(), "1H2M1H");

    }

    @Test ( enabled = false )
    public void testHardClipByReadCoordinates() {
        logger.warn("Executing testHardClipByReadCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReadCoordinates(0,3), new GATKSAMRecord(read.getHeader()));

        //clip 1 base at start
        expected = readClipper.hardClipByReadCoordinates(0,0);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,4));
        Assert.assertEquals(expected.getCigarString(), "1H3M");

        //clip 1 base at end
        expected = readClipper.hardClipByReadCoordinates(3,3);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,3).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,3));
        Assert.assertEquals(expected.getCigarString(), "3M1H");

        //clip 2 bases at start
        expected = readClipper.hardClipByReadCoordinates(0,1);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(2,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(2,4));
        Assert.assertEquals(expected.getCigarString(), "2H2M");

        //clip 2 bases at end
        expected = readClipper.hardClipByReadCoordinates(2,3);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,2).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,2));
        Assert.assertEquals(expected.getCigarString(), "2M2H");

    }

    @Test ( enabled = false )
    public void testHardClipByReferenceCoordinates() {
        logger.warn("Executing testHardClipByReferenceCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinates(1,4), new GATKSAMRecord(read.getHeader()));

        //clip 1 base at start
        expected = readClipper.hardClipByReferenceCoordinates(-1,1);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,4));
        Assert.assertEquals(expected.getCigarString(), "1H3M");

        //clip 1 base at end
        expected = readClipper.hardClipByReferenceCoordinates(3,-1);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,3).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,3));
        Assert.assertEquals(expected.getCigarString(), "3M1H");

        //clip 2 bases at start
        expected = readClipper.hardClipByReferenceCoordinates(-1,2);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(2,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(2,4));
        Assert.assertEquals(expected.getCigarString(), "2H2M");

        //clip 2 bases at end
        expected = readClipper.hardClipByReferenceCoordinates(2,-1);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,2).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,2));
        Assert.assertEquals(expected.getCigarString(), "2M2H");

    }

    @Test ( enabled = false )
    public void testHardClipByReferenceCoordinatesLeftTail() {
        logger.warn("Executing testHardClipByReferenceCoordinatesLeftTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesLeftTail(4), new GATKSAMRecord(read.getHeader()));

        //clip 1 base at start
        expected = readClipper.hardClipByReferenceCoordinatesLeftTail(1);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,4));
        Assert.assertEquals(expected.getCigarString(), "1H3M");

        //clip 2 bases at start
        expected = readClipper.hardClipByReferenceCoordinatesLeftTail(2);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(2,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(2,4));
        Assert.assertEquals(expected.getCigarString(), "2H2M");

    }

    @Test ( enabled = false )
    public void testHardClipByReferenceCoordinatesRightTail() {
        logger.warn("Executing testHardClipByReferenceCoordinatesRightTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesRightTail(1), new GATKSAMRecord(read.getHeader()));

        //clip 1 base at end
        expected = readClipper.hardClipByReferenceCoordinatesRightTail(3);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,3).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,3));
        Assert.assertEquals(expected.getCigarString(), "3M1H");

        //clip 2 bases at end
        expected = readClipper.hardClipByReferenceCoordinatesRightTail(2);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,2).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,2));
        Assert.assertEquals(expected.getCigarString(), "2M2H");

    }

    @Test ( enabled = false )
    public void testHardClipLowQualEnds() {
        logger.warn("Executing testHardClipByReferenceCoordinates");


        //Clip whole read
        Assert.assertEquals(readClipper.hardClipLowQualEnds((byte)64), new GATKSAMRecord(read.getHeader()));

        //clip 1 base at start
        expected = readClipper.hardClipLowQualEnds((byte)34);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,4));
        Assert.assertEquals(expected.getCigarString(), "1H3M");

        //clip 2 bases at start
        expected = readClipper.hardClipLowQualEnds((byte)44);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(2,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(2,4));
        Assert.assertEquals(expected.getCigarString(), "2H2M");

        // Reverse Quals sequence
        readClipper.getRead().setBaseQualityString("?5+!"); // 63,53,43,33

        //clip 1 base at end
        expected = readClipper.hardClipLowQualEnds((byte)34);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,3).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,3));
        Assert.assertEquals(expected.getCigarString(), "3M1H");

        //clip 2 bases at end
        expected = readClipper.hardClipLowQualEnds((byte)44);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(0,2).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(0,2));
        Assert.assertEquals(expected.getCigarString(), "2M2H");

        // revert Qual sequence
        readClipper.getRead().setBaseQualityString(QUALS);
    }
}
