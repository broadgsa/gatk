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
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 9/28/11
 * Time: 9:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadClipperUnitTest extends BaseTest {

    // TODO: Add error messages on failed tests


    //int debug = 0;

    GATKSAMRecord read, expected;
    ReadClipper readClipper;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?"; //ASCII values = 33,43,53,63

    // What the test read looks like
    // Ref:    1 2 3 4 5 6 7 8
    // Read:   0 1 2 3 - - - -
    // -----------------------------
    // Bases:  A C T G - - - -
    // Quals:  ! + 5 ? - - - -

    @BeforeTest
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, BASES.length());
        read.setReadBases(new String(BASES).getBytes());
        read.setBaseQualityString(new String(QUALS));

        readClipper = new ReadClipper(read);
        //logger.warn(read.getCigarString());
    }

    @Test ( enabled = true )
    public void testHardClipBothEndsByReferenceCoordinates() {
        init();
        logger.warn("Executing testHardClipBothEndsByReferenceCoordinates");
        //int debug = 1;
        //Clip whole read
        Assert.assertEquals(readClipper.hardClipBothEndsByReferenceCoordinates(1,1), new GATKSAMRecord(read.getHeader()));
        //clip 1 base
        expected = readClipper.hardClipBothEndsByReferenceCoordinates(1,4);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,3).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,3));
        Assert.assertEquals(expected.getCigarString(), "1H2M1H");

    }

    @Test ( enabled = false ) // TODO This fails at hardClipCigar and returns a NullPointerException
    public void testHardClipByReadCoordinates() {
        init();
        logger.warn("Executing testHardClipByReadCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReadCoordinates(0,3), new GATKSAMRecord(read.getHeader()));

        //clip 1 base at start
        System.out.println(readClipper.read.getCigarString());
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

    public void testIfEqual( GATKSAMRecord read, byte[] readBases, String baseQuals, String cigar) {
        Assert.assertEquals(read.getReadBases(), readBases);
        Assert.assertEquals(read.getBaseQualityString(), baseQuals);
        Assert.assertEquals(read.getCigarString(), cigar);
    }

    public class testParameter {
        int inputStart;
        int inputStop;
        int substringStart;
        int substringStop;
        String cigar;

        public testParameter(int InputStart, int InputStop, int SubstringStart, int SubstringStop, String Cigar) {
            inputStart = InputStart;
            inputStop = InputStop;
            substringStart = SubstringStart;
            substringStop = SubstringStop;
            cigar = Cigar;
            }
    }

    @Test ( enabled = true )
    public void testHardClipByReferenceCoordinates() {
        logger.warn("Executing testHardClipByReferenceCoordinates");
        //logger.warn(debug);
        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinates(1,4), new GATKSAMRecord(read.getHeader()));

        List<testParameter> testList = new LinkedList<testParameter>();
        testList.add(new testParameter(-1,1,1,4,"1H3M"));//clip 1 base at start
        testList.add(new testParameter(4,-1,0,3,"3M1H"));//clip 1 base at end
        testList.add(new testParameter(-1,2,2,4,"2H2M"));//clip 2 bases at start
        testList.add(new testParameter(3,-1,0,2,"2M2H"));//clip 2 bases at end

        for ( testParameter p : testList ) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.inputStop+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testIfEqual( readClipper.hardClipByReferenceCoordinates(p.inputStart,p.inputStop),
                    BASES.substring(p.substringStart,p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart,p.substringStop),
                    p.cigar );
        }

    }

    @Test ( enabled = true )
    public void testHardClipByReferenceCoordinatesLeftTail() {
        init();
        logger.warn("Executing testHardClipByReferenceCoordinatesLeftTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesLeftTail(4), new GATKSAMRecord(read.getHeader()));

        List<testParameter> testList = new LinkedList<testParameter>();
        testList.add(new testParameter(1,-1,1,4,"1H3M"));//clip 1 base at start
        testList.add(new testParameter(2,-1,2,4,"2H2M"));//clip 2 bases at start

        for ( testParameter p : testList ) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testIfEqual( readClipper.hardClipByReferenceCoordinatesLeftTail(p.inputStart),
                    BASES.substring(p.substringStart,p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart,p.substringStop),
                    p.cigar );
        }

    }

    @Test ( enabled = true )
    public void testHardClipByReferenceCoordinatesRightTail() {
        init();
        logger.warn("Executing testHardClipByReferenceCoordinatesRightTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesRightTail(1), new GATKSAMRecord(read.getHeader()));

        List<testParameter> testList = new LinkedList<testParameter>();
        testList.add(new testParameter(-1,4,0,3,"3M1H"));//clip 1 base at end
        testList.add(new testParameter(-1,3,0,2,"2M2H"));//clip 2 bases at end

        for ( testParameter p : testList ) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStop+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testIfEqual( readClipper.hardClipByReferenceCoordinatesRightTail(p.inputStop),
                    BASES.substring(p.substringStart,p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart,p.substringStop),
                    p.cigar );
        }

    }

    @Test ( enabled = false )  // TODO This function is returning null reads
    public void testHardClipLowQualEnds() {
        init();
        logger.warn("Executing testHardClipByReferenceCoordinates");


        //Clip whole read
        Assert.assertEquals(readClipper.hardClipLowQualEnds((byte)64), new GATKSAMRecord(read.getHeader()));

        //clip 1 base at start
        expected = readClipper.hardClipLowQualEnds((byte)34);
        logger.warn(expected.getBaseQualities().toString()+","+expected.getBaseQualityString());
        Assert.assertEquals(expected.getReadBases(), BASES.substring(1,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(1,4));
        Assert.assertEquals(expected.getCigarString(), "1H3M");

        //clip 2 bases at start
        expected = readClipper.hardClipLowQualEnds((byte)44);
        Assert.assertEquals(expected.getReadBases(), BASES.substring(2,4).getBytes());
        Assert.assertEquals(expected.getBaseQualityString(), QUALS.substring(2,4));
        Assert.assertEquals(expected.getCigarString(), "2H2M");

        // Reverse Quals sequence
        //readClipper.getRead().setBaseQualityString("?5+!"); // 63,53,43,33

        //clip 1 base at end
        expected = readClipper.hardClipLowQualEnds((byte)'!');
        logger.warn(expected.getBaseQualities().toString()+","+expected.getBaseQualityString());
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