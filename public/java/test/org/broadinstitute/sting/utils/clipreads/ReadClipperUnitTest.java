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

import net.sf.samtools.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 9/28/11
 * Time: 9:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadClipperUnitTest extends BaseTest {

    // TODO: exception testing, make cases that should fail will fail

    //int debug = 0;

    GATKSAMRecord read, expected;
    ReadClipper readClipper;
    final static String BASES = "ACTG";
    final static String QUALS = "!+5?"; //ASCII values = 33,43,53,63

    private void testBaseQualCigar(GATKSAMRecord read, byte[] readBases, byte[] baseQuals, String cigar) {
        // Because quals to char start at 33 for visibility
        baseQuals = subtractToArray(baseQuals, 33);

        Assert.assertEquals(read.getReadBases(), readBases);
        Assert.assertEquals(read.getBaseQualities(), baseQuals);
        Assert.assertEquals(read.getCigarString(), cigar);
    }

    private void testBaseQual(GATKSAMRecord read, byte[] readBases, byte[] baseQuals) {
        // Because quals to chars start at 33 for visibility
        baseQuals = subtractToArray(baseQuals, 33);

        if ( readBases.length > 0 && baseQuals.length > 0 ) {
            Assert.assertEquals(read.getReadBases(), readBases);
            Assert.assertEquals(read.getBaseQualities(), baseQuals);
        }
        else
            Assert.assertTrue(read.isEmpty());
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

    private byte[] subtractToArray (byte[] array, int n ) {
        if ( array == null)
            return null;

        byte[] output = new byte[array.length];

        for ( int i = 0; i < array.length; i++)
            output[i] = (byte) ( array[i] - n );

        return output;
    }

    // What the test read looks like
    // Ref:    1 2 3 4 5 6 7 8
    // Read:   0 1 2 3 - - - -
    // -----------------------------
    // Bases:  A C T G - - - -
    // Quals:  ! + 5 ? - - - -

    @BeforeMethod
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

    @Test ( enabled = true )
    public void testHardClipByReadCoordinates() {

        logger.warn("Executing testHardClipByReadCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReadCoordinates(0,3), new GATKSAMRecord(read.getHeader()));

        List<testParameter> testList = new LinkedList<testParameter>();
        testList.add(new testParameter(0,0,1,4,"1H3M"));//clip 1 base at start
        testList.add(new testParameter(3,3,0,3,"3M1H"));//clip 1 base at end
        testList.add(new testParameter(0,1,2,4,"2H2M"));//clip 2 bases at start
        testList.add(new testParameter(2,3,0,2,"2M2H"));//clip 2 bases at end

        for ( testParameter p : testList ) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.inputStop+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testBaseQualCigar(readClipper.hardClipByReadCoordinates(p.inputStart, p.inputStop),
                    BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
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
            testBaseQualCigar(readClipper.hardClipByReferenceCoordinates(p.inputStart, p.inputStop),
                    BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }

    }

    @Test ( enabled = true )
    public void testHardClipByReferenceCoordinatesLeftTail() {
        init();
        logger.warn("Executing testHardClipByReferenceCoordinatesLeftTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesLeftTail(4), new GATKSAMRecord(read.getHeader()));

        List<testParameter> testList = new LinkedList<testParameter>();
        testList.add(new testParameter(1, -1, 1, 4, "1H3M"));//clip 1 base at start
        testList.add(new testParameter(2, -1, 2, 4, "2H2M"));//clip 2 bases at start

        for ( testParameter p : testList ) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testBaseQualCigar(readClipper.hardClipByReferenceCoordinatesLeftTail(p.inputStart),
                    BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }

    }

    @Test ( enabled = true )
    public void testHardClipByReferenceCoordinatesRightTail() {
        init();
        logger.warn("Executing testHardClipByReferenceCoordinatesRightTail");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipByReferenceCoordinatesRightTail(1), new GATKSAMRecord(read.getHeader()));

        List<testParameter> testList = new LinkedList<testParameter>();
        testList.add(new testParameter(-1, 4, 0, 3, "3M1H"));//clip 1 base at end
        testList.add(new testParameter(-1, 3, 0, 2, "2M2H"));//clip 2 bases at end

        for ( testParameter p : testList ) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStop+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testBaseQualCigar(readClipper.hardClipByReferenceCoordinatesRightTail(p.inputStop),
                    BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart, p.substringStop).getBytes(),
                    p.cigar);
        }

    }

    @Test ( enabled = true )
    public void testHardClipLowQualEnds() {
        // Needs a thorough redesign
        logger.warn("Executing testHardClipByReferenceCoordinates");

        //Clip whole read
        Assert.assertEquals(readClipper.hardClipLowQualEnds((byte)64), new GATKSAMRecord(read.getHeader()));

        List<testParameter> testList = new LinkedList<testParameter>();
        testList.add(new testParameter(1,-1,1,4,"1H3M"));//clip 1 base at start
        testList.add(new testParameter(11,-1,2,4,"2H2M"));//clip 2 bases at start

        for ( testParameter p : testList ) {
            init();
            //logger.warn("Testing Parameters: " + p.inputStart+","+p.substringStart+","+p.substringStop+","+p.cigar);
            testBaseQualCigar(readClipper.hardClipLowQualEnds((byte) p.inputStart),
                    BASES.substring(p.substringStart, p.substringStop).getBytes(),
                    QUALS.substring(p.substringStart, p.substringStop).getBytes(),
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

  //  public class ReadMaker  {
    //Should have functions that can
    //    make basic read
    // make reads by cigar,
    // and by quals,
    // package in a readclipper
    // This has been already done in the current class
    //    GATKSAMRecord read;
    //    ReadClipper readClipper;

    public String cycleString ( String string, int length ) {
        String output = "";
        int cycles = ( length / string.length() ) + 1;

        for ( int i = 1; i < cycles; i++ )
            output += string;

        for ( int j = 0; output.length() < length; j++ )
            output += string.charAt( j % string.length() );

        return output;

    }

    private void makeFromCigar( Cigar cigar ) {
        readClipper = null;
        read = null;
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, cigar.getReadLength());
        read.setReadBases(cycleString(BASES, cigar.getReadLength()).getBytes());
        read.setBaseQualityString( cycleString(QUALS, cigar.getReadLength() ));
        read.setCigar(cigar);
        readClipper = new ReadClipper(read);
     }

 //   }

    private Set<Cigar> generateCigars() {

        // This function generates every permutation of cigar strings we need.

        LinkedHashSet<Cigar> output = new LinkedHashSet<Cigar>();

        List<Cigar> clippingOptionsStart = new LinkedList<Cigar>();
        clippingOptionsStart.add( new Cigar() );
        clippingOptionsStart.add( TextCigarCodec.getSingleton().decode("1H1S") );
        clippingOptionsStart.add( TextCigarCodec.getSingleton().decode("1S") );
        clippingOptionsStart.add( TextCigarCodec.getSingleton().decode("1H") );

        LinkedList<Cigar> clippingOptionsEnd = new LinkedList<Cigar>();
        clippingOptionsEnd.add( new Cigar() );
        clippingOptionsEnd.add( TextCigarCodec.getSingleton().decode("1S1H") );
        clippingOptionsEnd.add( TextCigarCodec.getSingleton().decode("1S"));
        clippingOptionsEnd.add( TextCigarCodec.getSingleton().decode("1H"));


        LinkedList<Cigar> indelOptions = new LinkedList<Cigar>();
        indelOptions.add( new Cigar() );
        //indelOptions.add( TextCigarCodec.getSingleton().decode("1I1D"));
        //indelOptions.add( TextCigarCodec.getSingleton().decode("1D1I") );
        indelOptions.add( TextCigarCodec.getSingleton().decode("1I") );
        indelOptions.add( TextCigarCodec.getSingleton().decode("1D") );

        // Start With M as base CigarElements, M,

        LinkedList<Cigar> base = new LinkedList<Cigar>();
        base.add( TextCigarCodec.getSingleton().decode("1M"));
        base.add( TextCigarCodec.getSingleton().decode("5M"));
        base.add( TextCigarCodec.getSingleton().decode("25M"));
        // Should indel be added as a base?

        // Nested loops W00t!
        for ( Cigar Base : base) {
            for ( Cigar indelStart: indelOptions) {
                for ( Cigar indelEnd: indelOptions) {
                    for ( Cigar clipStart: clippingOptionsStart) {
                        for ( Cigar clipEnd: clippingOptionsEnd) {
                            // Create a list of Cigar Elements and construct Cigar
                            List<CigarElement> CigarBuilder = new ArrayList<CigarElement>();
                            CigarBuilder.addAll(Base.getCigarElements());
                            CigarBuilder.addAll(indelStart.getCigarElements());
                            CigarBuilder.addAll(Base.getCigarElements());
                            CigarBuilder.addAll(indelEnd.getCigarElements());
                            CigarBuilder.addAll(clipEnd.getCigarElements());
                            //CigarBuilder.addAll(0, indelStart.getCigarElements());
                            CigarBuilder.addAll(0, clipStart.getCigarElements());
                            //System.out.println( new Cigar( removeConsecutiveElements(CigarBuilder) ).toString() );
                            output.add( new Cigar( removeConsecutiveElements(CigarBuilder) ) );

                        }
                    }
                }
            }
        }

        return output;
    }

    private List<CigarElement> removeConsecutiveElements(List<CigarElement> cigarBuilder) {
        LinkedList<CigarElement> output = new LinkedList<CigarElement>();
        for ( CigarElement E: cigarBuilder ) {
            if ( output.isEmpty() || output.getLast().getOperator() != E.getOperator())
                output.add(E);
        }
        return output;
    }

    @Test ( enabled = true )
    public void testHardClipSoftClippedBases() {

        // Generate a list of cigars to test
        for ( Cigar cigar: generateCigars() ) {
            //logger.warn("Testing Cigar: "+cigar.toString());
            makeFromCigar( cigar );


            try {
                readClipper.hardClipLeadingInsertions();
            }
            catch ( ReviewedStingException e ) {}


            int clipStart = 0;
            int clipEnd = 0;
            boolean expectEmptyRead = false;

            List<CigarElement> cigarElements = cigar.getCigarElements();
            int CigarListLength = cigarElements.size();


            // It will know what needs to be clipped based on the start and end of the string, hardclips and softclips
            // are added to the amount to clip
            if ( cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ) {
                //clipStart += cigarElements.get(0).getLength();
                if ( cigarElements.get(1).getOperator() == CigarOperator.SOFT_CLIP ) {
                    clipStart += cigarElements.get(1).getLength();
                    // Check for leading indel
                    if ( cigarElements.get(2).getOperator() == CigarOperator.INSERTION ) {
                        expectEmptyRead = true;
                    }
                }
                // Check for leading indel
                else if ( cigarElements.get(1).getOperator() == CigarOperator.INSERTION ) {
                    expectEmptyRead = true;
                }
            }
            else if ( cigarElements.get(0).getOperator() == CigarOperator.SOFT_CLIP ) {
                clipStart += cigarElements.get(0).getLength();
                // Check for leading indel
                if ( cigarElements.get(1).getOperator() == CigarOperator.INSERTION ) {
                    expectEmptyRead = true;
                }
            }
            //Check for leading indel
            else if ( cigarElements.get(0).getOperator() == CigarOperator.INSERTION ) {
                    expectEmptyRead = true;
                }

            if ( cigarElements.get(CigarListLength - 1).getOperator() == CigarOperator.HARD_CLIP ) {
                //clipEnd += cigarElements.get(CigarListLength - 1).getLength();
                if ( cigarElements.get(CigarListLength - 2).getOperator() == CigarOperator.SOFT_CLIP )
                    clipEnd += cigarElements.get(CigarListLength - 2).getLength();
            }
            else if ( cigarElements.get(CigarListLength - 1).getOperator() == CigarOperator.SOFT_CLIP )
                clipEnd += cigarElements.get(CigarListLength - 1).getLength();

            String readBases = readClipper.read.getReadString();
            String baseQuals = readClipper.read.getBaseQualityString();

            // "*" is the default empty-sequence-string and for our test it needs to be changed to  ""
            if (readBases.equals("*"))
                readBases = "";
            if (baseQuals.equals("*"))
                baseQuals = "";

            logger.warn(String.format("Testing cigar %s, expecting Base: %s and Qual: %s",
                    cigar.toString(), readBases.substring( clipStart, readBases.length() - clipEnd ),
                    baseQuals.substring( clipStart, baseQuals.length() - clipEnd ) ) );
            //if (expectEmptyRead)
            //    testBaseQual( readClipper.hardClipSoftClippedBases(), new byte[0], new byte[0] );
            //else
                testBaseQual(readClipper.hardClipSoftClippedBases(),
                    readBases.substring( clipStart, readBases.length() - clipEnd ).getBytes(),
                    baseQuals.substring( clipStart, baseQuals.length() - clipEnd ).getBytes());
            logger.warn("Cigar: "+cigar.toString()+" PASSED!");
        }
        // We will use testParameter in the following way
        // Right tail, left tail,

    }
}