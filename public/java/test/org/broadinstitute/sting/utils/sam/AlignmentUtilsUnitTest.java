/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.*;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class AlignmentUtilsUnitTest {
    private SAMFileHeader header;

    /** Basic aligned and mapped read. */
    private SAMRecord readMapped;

    /** Read with no contig specified in the read, -L UNMAPPED */
    private SAMRecord readNoReference;

    /** This read has a start position, but is flagged that it's not mapped. */
    private SAMRecord readUnmappedFlag;

    /** This read says it's aligned, but to a contig not in the header. */
    private SAMRecord readUnknownContig;

    /** This read says it's aligned, but actually has an unknown start. */
    private SAMRecord readUnknownStart;

    @BeforeClass
    public void init() {
        header = ArtificialSAMUtils.createArtificialSamHeader(3, 1, ArtificialSAMUtils.DEFAULT_READ_LENGTH * 2);

        readMapped = createMappedRead("mapped", 1);

        readNoReference = createUnmappedRead("unmappedNoReference");

        readUnmappedFlag = createMappedRead("unmappedFlagged", 2);
        readUnmappedFlag.setReadUnmappedFlag(true);

        readUnknownContig = createMappedRead("unknownContig", 3);
        readUnknownContig.setReferenceName("unknownContig");

        readUnknownStart = createMappedRead("unknownStart", 1);
        readUnknownStart.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
    }

    /**
     * Test for -L UNMAPPED
     */
    @DataProvider(name = "genomeLocUnmappedReadTests")
    public Object[][] getGenomeLocUnmappedReadTests() {
        return new Object[][] {
                new Object[] {readNoReference, true},
                new Object[] {readMapped, false},
                new Object[] {readUnmappedFlag, false},
                new Object[] {readUnknownContig, false},
                new Object[] {readUnknownStart, false}
        };
    }
    @Test(dataProvider = "genomeLocUnmappedReadTests")
    public void testIsReadGenomeLocUnmapped(SAMRecord read, boolean expected) {
        Assert.assertEquals(AlignmentUtils.isReadGenomeLocUnmapped(read), expected);
    }

    /**
     * Test for read being truly unmapped
     */
    @DataProvider(name = "unmappedReadTests")
    public Object[][] getUnmappedReadTests() {
        return new Object[][] {
                new Object[] {readNoReference, true},
                new Object[] {readMapped, false},
                new Object[] {readUnmappedFlag, true},
                new Object[] {readUnknownContig, false},
                new Object[] {readUnknownStart, true}
        };
    }
    @Test(dataProvider = "unmappedReadTests")
    public void testIsReadUnmapped(SAMRecord read, boolean expected) {
        Assert.assertEquals(AlignmentUtils.isReadUnmapped(read), expected);
    }

    private SAMRecord createUnmappedRead(String name) {
        return ArtificialSAMUtils.createArtificialRead(
                header,
                name,
                SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                SAMRecord.NO_ALIGNMENT_START,
                ArtificialSAMUtils.DEFAULT_READ_LENGTH);
    }

    private SAMRecord createMappedRead(String name, int start) {
        return ArtificialSAMUtils.createArtificialRead(
                header,
                name,
                0,
                start,
                ArtificialSAMUtils.DEFAULT_READ_LENGTH);
    }

    private final List<List<CigarElement>> makeCigarElementCombinations() {
        // this functionality can be adapted to provide input data for whatever you might want in your data
        final List<CigarElement> cigarElements = new LinkedList<CigarElement>();
        for ( final int size : Arrays.asList(0, 10) ) {
            for ( final CigarOperator op : CigarOperator.values() ) {
                cigarElements.add(new CigarElement(size, op));
            }
        }

        final List<List<CigarElement>> combinations = new LinkedList<List<CigarElement>>();
        for ( final int nElements : Arrays.asList(1, 2, 3) ) {
            combinations.addAll(Utils.makePermutations(cigarElements, nElements, true));
        }

        return combinations;
    }


    @DataProvider(name = "NumAlignedBasesCountingSoftClips")
    public Object[][] makeNumAlignedBasesCountingSoftClips() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final EnumSet<CigarOperator> alignedToGenome = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X, CigarOperator.S);
        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            int n = 0;
            for ( final CigarElement elt : elements ) n += alignedToGenome.contains(elt.getOperator()) ? elt.getLength() : 0;
            tests.add(new Object[]{new Cigar(elements), n});
        }

        tests.add(new Object[]{null, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "NumAlignedBasesCountingSoftClips")
    public void testNumAlignedBasesCountingSoftClips(final Cigar cigar, final int expected) {
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, cigar == null ? 10 : cigar.getReadLength());
        read.setCigar(cigar);
        Assert.assertEquals(AlignmentUtils.getNumAlignedBasesCountingSoftClips(read), expected, "Cigar " + cigar + " failed NumAlignedBasesCountingSoftClips");
    }

    @DataProvider(name = "CigarHasZeroElement")
    public Object[][] makeCigarHasZeroElement() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            boolean hasZero = false;
            for ( final CigarElement elt : elements ) hasZero = hasZero || elt.getLength() == 0;
            tests.add(new Object[]{new Cigar(elements), hasZero});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CigarHasZeroElement")
    public void testCigarHasZeroSize(final Cigar cigar, final boolean hasZero) {
        Assert.assertEquals(AlignmentUtils.cigarHasZeroSizeElement(cigar), hasZero, "Cigar " + cigar.toString() + " failed cigarHasZeroSizeElement");
    }

    @DataProvider(name = "NumHardClipped")
    public Object[][] makeNumHardClipped() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            int n = 0;
            for ( final CigarElement elt : elements ) n += elt.getOperator() == CigarOperator.H ? elt.getLength() : 0;
            tests.add(new Object[]{new Cigar(elements), n});
        }

        tests.add(new Object[]{null, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "NumHardClipped")
    public void testNumHardClipped(final Cigar cigar, final int expected) {
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, cigar == null ? 10 : cigar.getReadLength());
        read.setCigar(cigar);
        Assert.assertEquals(AlignmentUtils.getNumHardClippedBases(read), expected, "Cigar " + cigar + " failed num hard clips");
    }

    @DataProvider(name = "NumAlignedBlocks")
    public Object[][] makeNumAlignedBlocks() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            int n = 0;
            for ( final CigarElement elt : elements ) {
                switch ( elt.getOperator() ) {
                    case M:case X:case EQ: n++; break;
                    default: break;
                }
            }
            tests.add(new Object[]{new Cigar(elements), n});
        }

        tests.add(new Object[]{null, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "NumAlignedBlocks")
    public void testNumAlignedBlocks(final Cigar cigar, final int expected) {
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, cigar == null ? 10 : cigar.getReadLength());
        read.setCigar(cigar);
        Assert.assertEquals(AlignmentUtils.getNumAlignmentBlocks(read), expected, "Cigar " + cigar + " failed NumAlignedBlocks");
    }

    @Test
    public void testConsolidateCigar() {
        {
            //1M1M1M1D2M1M --> 3M1D3M
            List<CigarElement> list = new ArrayList<CigarElement>();
            list.add( new CigarElement(1, CigarOperator.M));
            list.add( new CigarElement(1, CigarOperator.M));
            list.add( new CigarElement(1, CigarOperator.M));
            list.add( new CigarElement(1, CigarOperator.D));
            list.add( new CigarElement(2, CigarOperator.M));
            list.add( new CigarElement(1, CigarOperator.M));
            Cigar unconsolidatedCigar = new Cigar(list);

            list.clear();
            list.add( new CigarElement(3, CigarOperator.M));
            list.add( new CigarElement(1, CigarOperator.D));
            list.add( new CigarElement(3, CigarOperator.M));
            Cigar consolidatedCigar = new Cigar(list);

            Assert.assertEquals(consolidatedCigar.toString(), AlignmentUtils.consolidateCigar(unconsolidatedCigar).toString());
        }

        {
            //6M6M6M --> 18M
            List<CigarElement> list = new ArrayList<CigarElement>();
            list.add( new CigarElement(6, CigarOperator.M));
            list.add( new CigarElement(6, CigarOperator.M));
            list.add( new CigarElement(6, CigarOperator.M));
            Cigar unconsolidatedCigar = new Cigar(list);

            list.clear();
            list.add( new CigarElement(18, CigarOperator.M));
            Cigar consolidatedCigar = new Cigar(list);

            Assert.assertEquals(consolidatedCigar.toString(), AlignmentUtils.consolidateCigar(unconsolidatedCigar).toString());
        }
    }

    @DataProvider(name = "SoftClipsDataProvider")
    public Object[][] makeSoftClipsDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final int lengthOfLeftClip : Arrays.asList(0, 1, 10) ) {
            for ( final int lengthOfRightClip : Arrays.asList(0, 1, 10) ) {
                for ( final int qualThres : Arrays.asList(10, 20, 30) ) {
                    for ( final String middleOp : Arrays.asList("M", "D") ) {
                        for ( final int matchSize : Arrays.asList(0, 1, 10) ) {
                            final byte[] left = makeQualArray(lengthOfLeftClip, qualThres);
                            final byte[] right = makeQualArray(lengthOfRightClip, qualThres);
                            int n = 0;
                            for ( int i = 0; i < left.length; i++ ) n += left[i] > qualThres ? 1 : 0;
                            for ( int i = 0; i < right.length; i++ ) n += right[i] > qualThres ? 1 : 0;
                            tests.add(new Object[]{left, matchSize, middleOp, right, qualThres, n});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private byte[] makeQualArray(final int length, final int qualThreshold) {
        final byte[] array = new byte[length];
        for ( int i = 0; i < array.length; i++ )
            array[i] = (byte)(qualThreshold + ( i % 2 == 0 ? 1 : - 1 ));
        return array;
    }

    @Test(dataProvider = "SoftClipsDataProvider")
    public void testSoftClipsData(final byte[] qualsOfSoftClipsOnLeft, final int middleSize, final String middleOp, final byte[] qualOfSoftClipsOnRight, final int qualThreshold, final int numExpected) {
        final int readLength = (middleOp.equals("D") ? 0 : middleSize) + qualOfSoftClipsOnRight.length + qualsOfSoftClipsOnLeft.length;
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
        final byte[] bases = Utils.dupBytes((byte) 'A', readLength);
        final byte[] matchBytes = middleOp.equals("D") ? new byte[]{} : Utils.dupBytes((byte)30, middleSize);
        final byte[] quals = ArrayUtils.addAll(ArrayUtils.addAll(qualsOfSoftClipsOnLeft, matchBytes), qualOfSoftClipsOnRight);

        // set the read's bases and quals
        read.setReadBases(bases);
        read.setBaseQualities(quals);

        final StringBuilder cigar = new StringBuilder();
        if (qualsOfSoftClipsOnLeft.length > 0 ) cigar.append(qualsOfSoftClipsOnLeft.length + "S");
        if (middleSize > 0 ) cigar.append(middleSize + middleOp);
        if (qualOfSoftClipsOnRight.length > 0 ) cigar.append(qualOfSoftClipsOnRight.length + "S");

        read.setCigarString(cigar.toString());

        final int actual = AlignmentUtils.calcNumHighQualitySoftClips(read, (byte) qualThreshold);
        Assert.assertEquals(actual, numExpected, "Wrong number of soft clips detected for read " + read.getSAMString());
    }
}
