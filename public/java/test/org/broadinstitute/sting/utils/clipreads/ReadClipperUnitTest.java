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
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
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
        logger.warn("Executing testHardClipLowQualEnds");

        // Testing clipping that ends inside an insertion
        final byte[] BASES = {'A','C','G','T','A','C','G','T'};
        final byte[] QUALS = {2, 2, 2, 2, 20, 20, 20, 2};
        final String CIGAR = "1S1M5I1S";

        final byte[] CLIPPED_BASES = {};
        final byte[] CLIPPED_QUALS = {};
        final String CLIPPED_CIGAR = "";


        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(BASES, QUALS, CIGAR);
        GATKSAMRecord expected = ArtificialSAMUtils.createArtificialRead(CLIPPED_BASES, CLIPPED_QUALS, CLIPPED_CIGAR);

        ReadClipper lowQualClipper = new ReadClipper(read);
        ClipReadsTestUtils.assertEqualReads(lowQualClipper.hardClipLowQualEnds((byte) 2), expected);


    }

    @Test(enabled = true)
    public void testHardClipSoftClippedBases() {

        // Generate a list of cigars to test
        for (Cigar cigar : ClipReadsTestUtils.generateCigars()) {
            GATKSAMRecord read = ClipReadsTestUtils.makeReadFromCigar(cigar);
            readClipper = new ReadClipper(read);
            GATKSAMRecord clippedRead = readClipper.hardClipSoftClippedBases();

            int sumHardClips = 0;
            int sumMatches = 0;

            boolean tail = true;
            for (CigarElement element : read.getCigar().getCigarElements()) {
                // Assuming cigars are well formed, if we see S or H, it means we're on the tail (left or right)
                if (element.getOperator() == CigarOperator.HARD_CLIP || element.getOperator() == CigarOperator.SOFT_CLIP)
                    tail = true;

                // Adds all H, S and D's (next to hard/soft clips).
                // All these should be hard clips after clipping.
                if (tail && (element.getOperator() == CigarOperator.HARD_CLIP || element.getOperator() == CigarOperator.SOFT_CLIP || element.getOperator() == CigarOperator.DELETION))
                    sumHardClips += element.getLength();

                // this means we're no longer on the tail (insertions can still potentially be the tail because
                // of the current contract of clipping out hanging insertions
                else if (element.getOperator() != CigarOperator.INSERTION)
                    tail = false;

                // Adds all matches to verify that they remain the same after clipping
                if (element.getOperator() == CigarOperator.MATCH_OR_MISMATCH)
                    sumMatches += element.getLength();
            }

            for (CigarElement element : clippedRead.getCigar().getCigarElements()) {
                // Test if clipped read has Soft Clips (shouldn't have any!)
                Assert.assertTrue( element.getOperator() != CigarOperator.SOFT_CLIP, String.format("Cigar %s -> %s -- FAILED (resulting cigar has soft clips)", read.getCigarString(), clippedRead.getCigarString()));

                // Keep track of the total number of Hard Clips after clipping to make sure everything was accounted for
                if (element.getOperator() == CigarOperator.HARD_CLIP)
                    sumHardClips -= element.getLength();

                // Make sure all matches are still there
                if (element.getOperator() == CigarOperator.MATCH_OR_MISMATCH)
                    sumMatches -= element.getLength();
            }
            Assert.assertTrue( sumHardClips == 0, String.format("Cigar %s -> %s -- FAILED (number of hard clips mismatched by %d)", read.getCigarString(), clippedRead.getCigarString(), sumHardClips));
            Assert.assertTrue( sumMatches == 0, String.format("Cigar %s -> %s -- FAILED (number of matches mismatched by %d)", read.getCigarString(), clippedRead.getCigarString(), sumMatches));


            logger.warn(String.format("Cigar %s -> %s -- PASSED!", read.getCigarString(), clippedRead.getCigarString()));
        }
    }
}