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

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.TextCigarCodec;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
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

    // TODO: exception testing, make cases that should fail will fail

    // TODO: add indels to all test cases

    ReadClipper readClipper;

    @BeforeMethod
    public void init() {
        readClipper = new ReadClipper(ClipReadsTestUtils.makeRead());
    }

    @Test(enabled = true)
    public void testHardClipBothEndsByReferenceCoordinates() {

        logger.warn("Executing testHardClipBothEndsByReferenceCoordinates");
        //int debug = 1;
        //Clip whole read
        Assert.assertEquals(readClipper.hardClipBothEndsByReferenceCoordinates(10, 10), new GATKSAMRecord(readClipper.read.getHeader()));

        //clip 1 base
        ClipReadsTestUtils.testBaseQualCigar(readClipper.hardClipBothEndsByReferenceCoordinates(10, 13),
                ClipReadsTestUtils.BASES.substring(1, 3).getBytes(), ClipReadsTestUtils.QUALS.substring(1, 3).getBytes(),
                "1H2M1H");

        List<CigarStringTestPair> cigarStringTestPairs = new LinkedList<CigarStringTestPair>();
        cigarStringTestPairs.add(new CigarStringTestPair("5M1D1M2I4M5I6M1D3M2I100M", "1H4M1D1M2I4M5I6M1D3M2I99M1H"));
        //cigarStringTestPairs.add( new CigarStringTestPair("5M1I1M2I1M","1H4M1I1M2I1H"));
        cigarStringTestPairs.add(new CigarStringTestPair("1S1M1I1M1I1M1I1M1I1M1I1M1S", "1H1M1I1M1I1M1I1M1I1M1I1M1H"));
        cigarStringTestPairs.add(new CigarStringTestPair("1S1M1D1M1D1M1D1M1D1M1D1M1S", "1H1M1D1M1D1M1D1M1D1M1D1M1H"));

        for (CigarStringTestPair pair : cigarStringTestPairs) {
            readClipper = new ReadClipper(ClipReadsTestUtils.makeReadFromCigar(TextCigarCodec.getSingleton().decode(pair.toTest)));
            ClipReadsTestUtils.testCigar(readClipper.hardClipBothEndsByReferenceCoordinates(
                    ReadUtils.getRefCoordSoftUnclippedStart(readClipper.read),
                    ReadUtils.getRefCoordSoftUnclippedEnd(readClipper.read)),
                    pair.expected);
        }
        /*
        for ( Cigar cigar: ClipReadsTestUtils.generateCigars() ) {
            // The read has to be long enough to clip one base from each side
            // This also filters a lot of cigars
            if ( cigar.getReadLength() > 26 ) {
                readClipper = new ReadClipper(ClipReadsTestUtils.makeReadFromCigar( cigar ));
                System.out.println( "Testing Cigar: "+cigar.toString() ) ;
                //cigar length reference plus soft clip

                ClipReadsTestUtils.testBaseQual(
                        readClipper.hardClipBothEndsByReferenceCoordinates(
                                                              ReadUtils.getRefCoordSoftUnclippedStart(readClipper.read),
                                                              ReadUtils.getRefCoordSoftUnclippedEnd(readClipper.read) ),
                        readClipper.read.getReadString().substring(1, (cigar.getReadLength() - 1)).getBytes(),
                        readClipper.read.getBaseQualityString().substring(1, (cigar.getReadLength() - 1)).getBytes());
            }
        }
        */
    }

    @Test(enabled = true)
    public void testHardClipByReadCoordinates() {

        logger.warn("Executing testHardClipByReadCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReadCoordinates(0, 3), new GATKSAMRecord(readClipper.read.getHeader()));

        List<TestParameter> testList = new LinkedList<TestParameter>();
        testList.add(new TestParameter(0, 0, 1, 4, "1H3M"));//clip 1 base at start
        testList.add(new TestParameter(3, 3, 0, 3, "3M1H"));//clip 1 base at end
        testList.add(new TestParameter(0, 1, 2, 4, "2H2M"));//clip 2 bases at start
        testList.add(new TestParameter(2, 3, 0, 2, "2M2H"));//clip 2 bases at end

        for (TestParameter p : testList) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.inputStop+","+p.substringStart+","+p.substringStop+","+p.cigar);
            ClipReadsTestUtils.testBaseQualCigar(readClipper.hardClipByReadCoordinates(p.inputStart, p.inputStop),
                    ClipReadsTestUtils.BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    ClipReadsTestUtils.QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }

    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinates() {
        logger.warn("Executing testHardClipByReferenceCoordinates");
        //logger.warn(debug);
        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinates(10, 13), new GATKSAMRecord(readClipper.read.getHeader()));

        List<TestParameter> testList = new LinkedList<TestParameter>();
        testList.add(new TestParameter(-1, 10, 1, 4, "1H3M"));//clip 1 base at start
        testList.add(new TestParameter(13, -1, 0, 3, "3M1H"));//clip 1 base at end
        testList.add(new TestParameter(-1, 11, 2, 4, "2H2M"));//clip 2 bases at start
        testList.add(new TestParameter(12, -1, 0, 2, "2M2H"));//clip 2 bases at end

        for (TestParameter p : testList) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.inputStop+","+p.substringStart+","+p.substringStop+","+p.cigar);
            ClipReadsTestUtils.testBaseQualCigar(readClipper.hardClipByReferenceCoordinates(p.inputStart, p.inputStop),
                    ClipReadsTestUtils.BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    ClipReadsTestUtils.QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }

        List<CigarStringTestPair> cigarStringTestPairs = new LinkedList<CigarStringTestPair>();
        cigarStringTestPairs.add(new CigarStringTestPair("5M1D1M2I4M5I6M1D3M2I100M", "1H4M1D1M2I4M5I6M1D3M2I100M"));
        //cigarStringTestPairs.add( new CigarStringTestPair("5M1I1M2I1M","1H4M1I1M2I1M"));
        cigarStringTestPairs.add(new CigarStringTestPair("1S1M1I1M1I1M1I1M1I1M1I1M1S", "1H1M1I1M1I1M1I1M1I1M1I1M1S"));
        cigarStringTestPairs.add(new CigarStringTestPair("1S1M1D1M1D1M1D1M1D1M1D1M1S", "1H1M1D1M1D1M1D1M1D1M1D1M1S"));

        //Clips only first base
        for (CigarStringTestPair pair : cigarStringTestPairs) {
            readClipper = new ReadClipper(ClipReadsTestUtils.makeReadFromCigar(TextCigarCodec.getSingleton().decode(pair.toTest)));
            ClipReadsTestUtils.testCigar(readClipper.hardClipByReadCoordinates(0, 0), pair.expected);
        }
        /*
        for ( Cigar cigar: ClipReadsTestUtils.generateCigars() ) {
            // The read has to be long enough to clip one base
            // This also filters a lot of cigars
            if ( cigar.getReadLength() > 26 ) {
                readClipper = new ReadClipper(ClipReadsTestUtils.makeReadFromCigar( cigar ));
                System.out.println( "Testing Cigar: "+cigar.toString() ) ;
                //cigar length reference plus soft clip

                // Clip first read
                ClipReadsTestUtils.testBaseQual(
                        readClipper.hardClipByReadCoordinates(0,0),
                        readClipper.read.getReadString().substring(1, cigar.getReadLength()).getBytes(),
                        readClipper.read.getBaseQualityString().substring(1, cigar.getReadLength()).getBytes());
            }
        }
        */
    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinatesLeftTail() {
        init();
        logger.warn("Executing testHardClipByReferenceCoordinatesLeftTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesLeftTail(13), new GATKSAMRecord(readClipper.read.getHeader()));

        List<TestParameter> testList = new LinkedList<TestParameter>();
        testList.add(new TestParameter(10, -1, 1, 4, "1H3M"));//clip 1 base at start
        testList.add(new TestParameter(11, -1, 2, 4, "2H2M"));//clip 2 bases at start

        for (TestParameter p : testList) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.substringStart+","+p.substringStop+","+p.cigar);
            ClipReadsTestUtils.testBaseQualCigar(readClipper.hardClipByReferenceCoordinatesLeftTail(p.inputStart),
                    ClipReadsTestUtils.BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    ClipReadsTestUtils.QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }

    }

    @Test(enabled = true)
    public void testHardClipByReferenceCoordinatesRightTail() {
        init();
        logger.warn("Executing testHardClipByReferenceCoordinatesRightTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesRightTail(10), new GATKSAMRecord(readClipper.read.getHeader()));

        List<TestParameter> testList = new LinkedList<TestParameter>();
        testList.add(new TestParameter(-1, 13, 0, 3, "3M1H"));//clip 1 base at end
        testList.add(new TestParameter(-1, 12, 0, 2, "2M2H"));//clip 2 bases at end

        for (TestParameter p : testList) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStop+","+p.substringStart+","+p.substringStop+","+p.cigar);
            ClipReadsTestUtils.testBaseQualCigar(readClipper.hardClipByReferenceCoordinatesRightTail(p.inputStop),
                    ClipReadsTestUtils.BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    ClipReadsTestUtils.QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }

    }

    @Test(enabled = true)
    public void testHardClipLowQualEnds() {
        // Needs a thorough redesign
        logger.warn("Executing testHardClipByReferenceCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipLowQualEnds((byte) 64), new GATKSAMRecord(readClipper.read.getHeader()));

        List<TestParameter> testList = new LinkedList<TestParameter>();
        testList.add(new TestParameter(1, -1, 1, 4, "1H3M"));//clip 1 base at start
        testList.add(new TestParameter(11, -1, 2, 4, "2H2M"));//clip 2 bases at start

        for (TestParameter p : testList) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.substringStart+","+p.substringStop+","+p.cigar);
            ClipReadsTestUtils.testBaseQualCigar(readClipper.hardClipLowQualEnds((byte) p.inputStart),
                    ClipReadsTestUtils.BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    ClipReadsTestUtils.QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }
        /*      todo find a better way to test lowqual tail clipping on both sides
        // Reverse Quals sequence
        readClipper.getRead().setBaseQualityString("?5+!"); // 63,53,43,33

        testList = new LinkedList<testParameter>();
        testList.add(new testParameter(1,-1,0,3,"3M1H"));//clip 1 base at end
        testList.add(new testParameter(11,-1,0,2,"2M2H"));//clip 2 bases at end

        for ( testParameter p : testList ) {
            init();
            readClipper.getRead().setBaseQualityString("?5+!"); // 63,53,43,33
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testBaseQualCigar( readClipper.hardClipLowQualEnds( (byte)p.inputStart ),
                    BASES.substring(p.substringStart,p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart,p.substringStop),
                    p.cigar );
        }
        */
    }

    @Test(enabled = false)
    public void testHardClipSoftClippedBases() {

        // Generate a list of cigars to test
        for (Cigar cigar : ClipReadsTestUtils.generateCigars()) {
            //logger.warn("Testing Cigar: "+cigar.toString());
            readClipper = new ReadClipper(ClipReadsTestUtils.makeReadFromCigar(cigar));

            int clipStart = 0;
            int clipEnd = 0;
            boolean expectEmptyRead = false;

            List<CigarElement> cigarElements = cigar.getCigarElements();
            int CigarListLength = cigarElements.size();

            // It will know what needs to be clipped based on the start and end of the string, hardclips and softclips
            // are added to the amount to clip
            if (cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP) {
                //clipStart += cigarElements.get(0).getLength();
                if (cigarElements.get(1).getOperator() == CigarOperator.SOFT_CLIP) {
                    clipStart += cigarElements.get(1).getLength();
                    // Check for leading indel
                    if (cigarElements.get(2).getOperator() == CigarOperator.INSERTION) {
                        expectEmptyRead = true;
                    }
                }
                // Check for leading indel
                else if (cigarElements.get(1).getOperator() == CigarOperator.INSERTION) {
                    expectEmptyRead = true;
                }
            } else if (cigarElements.get(0).getOperator() == CigarOperator.SOFT_CLIP) {
                clipStart += cigarElements.get(0).getLength();
                // Check for leading indel
                if (cigarElements.get(1).getOperator() == CigarOperator.INSERTION) {
                    expectEmptyRead = true;
                }
            }
            //Check for leading indel
            else if (cigarElements.get(0).getOperator() == CigarOperator.INSERTION) {
                expectEmptyRead = true;
            }

            if (cigarElements.get(CigarListLength - 1).getOperator() == CigarOperator.HARD_CLIP) {
                //clipEnd += cigarElements.get(CigarListLength - 1).getLength();
                if (cigarElements.get(CigarListLength - 2).getOperator() == CigarOperator.SOFT_CLIP)
                    clipEnd += cigarElements.get(CigarListLength - 2).getLength();
            } else if (cigarElements.get(CigarListLength - 1).getOperator() == CigarOperator.SOFT_CLIP)
                clipEnd += cigarElements.get(CigarListLength - 1).getLength();

            String readBases = readClipper.read.getReadString();
            String baseQuals = readClipper.read.getBaseQualityString();

            // "*" is the default empty-sequence-string and for our test it needs to be changed to  ""
            if (readBases.equals("*"))
                readBases = "";
            if (baseQuals.equals("*"))
                baseQuals = "";

            logger.warn(String.format("Testing cigar %s, expecting Base: %s and Qual: %s",
                    cigar.toString(), readBases.substring(clipStart, readBases.length() - clipEnd),
                    baseQuals.substring(clipStart, baseQuals.length() - clipEnd)));
            //if (expectEmptyRead)
            //    testBaseQual( readClipper.hardClipSoftClippedBases(), new byte[0], new byte[0] );
            //else
            ClipReadsTestUtils.testBaseQual(readClipper.hardClipSoftClippedBases(),
                    readBases.substring(clipStart, readBases.length() - clipEnd).getBytes(),
                    baseQuals.substring(clipStart, baseQuals.length() - clipEnd).getBytes());
            logger.warn("Cigar: " + cigar.toString() + " PASSED!");
        }
        // We will use testParameter in the following way
        // Right tail, left tail,

    }
}