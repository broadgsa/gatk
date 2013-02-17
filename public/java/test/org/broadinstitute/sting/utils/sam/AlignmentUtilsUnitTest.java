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
import org.broadinstitute.sting.utils.pileup.PileupElement;
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

    ////////////////////////////////////////////
    // Test AlignmentUtils.getMismatchCount() //
    ////////////////////////////////////////////

    @DataProvider(name = "MismatchCountDataProvider")
    public Object[][] makeMismatchCountDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final int readLength = 20;
        final int lengthOfIndel = 2;
        final int locationOnReference = 10;
        final byte[] reference = Utils.dupBytes((byte)'A', readLength);
        final byte[] quals = Utils.dupBytes((byte)'A', readLength);


        for ( int startOnRead = 0; startOnRead <= readLength; startOnRead++ ) {
            for ( int basesToRead = 0; basesToRead <= readLength; basesToRead++ ) {
                for ( final int lengthOfSoftClip : Arrays.asList(0, 1, 10) ) {
                    for ( final int lengthOfFirstM : Arrays.asList(0, 3) ) {
                        for ( final char middleOp : Arrays.asList('M', 'D', 'I') ) {
                            for ( final int mismatchLocation : Arrays.asList(-1, 0, 5, 10, 15, 19) ) {

                                final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, locationOnReference, readLength);

                                // set the read's bases and quals
                                final byte[] readBases = reference.clone();
                                // create the mismatch if requested
                                if ( mismatchLocation != -1 )
                                    readBases[mismatchLocation] = (byte)'C';
                                read.setReadBases(readBases);
                                read.setBaseQualities(quals);

                                // create the CIGAR string
                                read.setCigarString(buildTestCigarString(middleOp, lengthOfSoftClip, lengthOfFirstM, lengthOfIndel, readLength));

                                // now, determine whether or not there's a mismatch
                                final boolean isMismatch;
                                if ( mismatchLocation < startOnRead || mismatchLocation >= startOnRead + basesToRead || mismatchLocation < lengthOfSoftClip ) {
                                    isMismatch = false;
                                } else if ( middleOp == 'M' || middleOp == 'D' || mismatchLocation < lengthOfSoftClip + lengthOfFirstM || mismatchLocation >= lengthOfSoftClip + lengthOfFirstM + lengthOfIndel ) {
                                    isMismatch = true;
                                } else {
                                    isMismatch = false;
                                }

                                tests.add(new Object[]{read, locationOnReference, startOnRead, basesToRead, isMismatch});
                            }
                        }
                    }
                }
            }
        }

        // Adding test to make sure soft-clipped reads go through the exceptions thrown at the beginning of the getMismatchCount method
        // todo: incorporate cigars with right-tail soft-clips in the systematic tests above.
        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 10, 20);
        read.setReadBases(reference);
        read.setBaseQualities(quals);
        read.setCigarString("10S5M5S");
        tests.add(new Object[]{read, 10, read.getAlignmentStart(), read.getReadLength(), false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MismatchCountDataProvider")
    public void testMismatchCountData(final GATKSAMRecord read, final int refIndex, final int startOnRead, final int basesToRead, final boolean isMismatch) {
        final byte[] reference = Utils.dupBytes((byte)'A', 100);
        final int actual = AlignmentUtils.getMismatchCount(read, reference, refIndex, startOnRead, basesToRead).numMismatches;
        Assert.assertEquals(actual, isMismatch ? 1 : 0, "Wrong number of mismatches detected for read " + read.getSAMString());
    }

    private static String buildTestCigarString(final char middleOp, final int lengthOfSoftClip, final int lengthOfFirstM, final int lengthOfIndel, final int readLength) {
        final StringBuilder cigar = new StringBuilder();
        int remainingLength = readLength;

        // add soft clips to the beginning of the read
        if (lengthOfSoftClip > 0 ) {
            cigar.append(lengthOfSoftClip).append("S");
            remainingLength -= lengthOfSoftClip;
        }

        if ( middleOp == 'M' ) {
            cigar.append(remainingLength).append("M");
        } else {
            if ( lengthOfFirstM > 0 ) {
                cigar.append(lengthOfFirstM).append("M");
                remainingLength -= lengthOfFirstM;
            }

            if ( middleOp == 'D' ) {
                cigar.append(lengthOfIndel).append("D");
            } else {
                cigar.append(lengthOfIndel).append("I");
                remainingLength -= lengthOfIndel;
            }
            cigar.append(remainingLength).append("M");
        }

        return cigar.toString();
    }

    ////////////////////////////////////////////////////////
    // Test AlignmentUtils.calcAlignmentByteArrayOffset() //
    ////////////////////////////////////////////////////////

    @DataProvider(name = "AlignmentByteArrayOffsetDataProvider")
    public Object[][] makeAlignmentByteArrayOffsetDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final int readLength = 20;
        final int lengthOfIndel = 2;
        final int locationOnReference = 20;

        for ( int offset = 0; offset < readLength; offset++ ) {
            for ( final int lengthOfSoftClip : Arrays.asList(0, 1, 10) ) {
                for ( final int lengthOfFirstM : Arrays.asList(0, 3) ) {
                    for ( final char middleOp : Arrays.asList('M', 'D', 'I') ) {

                        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, locationOnReference, readLength);
                        // create the CIGAR string
                        read.setCigarString(buildTestCigarString(middleOp, lengthOfSoftClip, lengthOfFirstM, lengthOfIndel, readLength));

                        // now, determine the expected alignment offset
                        final int expected;
                        boolean isDeletion = false;
                        if ( offset < lengthOfSoftClip ) {
                            expected = 0;
                        } else if ( middleOp == 'M' || offset < lengthOfSoftClip + lengthOfFirstM ) {
                            expected = offset - lengthOfSoftClip;
                        } else if ( offset < lengthOfSoftClip + lengthOfFirstM + lengthOfIndel ) {
                            if ( middleOp == 'D' ) {
                                isDeletion = true;
                                expected = offset - lengthOfSoftClip;
                            } else {
                                expected = lengthOfFirstM;
                            }
                        } else {
                            expected = offset - lengthOfSoftClip - (middleOp == 'I' ? lengthOfIndel : -lengthOfIndel);
                        }

                        tests.add(new Object[]{read.getCigar(), offset, expected, isDeletion, lengthOfSoftClip});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AlignmentByteArrayOffsetDataProvider")
    public void testAlignmentByteArrayOffsetData(final Cigar cigar, final int offset, final int expectedResult, final boolean isDeletion, final int lengthOfSoftClip) {
        final int actual = AlignmentUtils.calcAlignmentByteArrayOffset(cigar, isDeletion ? -1 : offset, isDeletion, 20, 20 + offset - lengthOfSoftClip);
        Assert.assertEquals(actual, expectedResult, "Wrong alignment offset detected for cigar " + cigar.toString());
    }

    ////////////////////////////////////////////////////
    // Test AlignmentUtils.readToAlignmentByteArray() //
    ////////////////////////////////////////////////////

    @DataProvider(name = "ReadToAlignmentByteArrayDataProvider")
    public Object[][] makeReadToAlignmentByteArrayDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final int readLength = 20;
        final int lengthOfIndel = 2;
        final int locationOnReference = 20;

        for ( final int lengthOfSoftClip : Arrays.asList(0, 1, 10) ) {
            for ( final int lengthOfFirstM : Arrays.asList(0, 3) ) {
                for ( final char middleOp : Arrays.asList('M', 'D', 'I') ) {

                    final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, locationOnReference, readLength);
                    // create the CIGAR string
                    read.setCigarString(buildTestCigarString(middleOp, lengthOfSoftClip, lengthOfFirstM, lengthOfIndel, readLength));

                    // now, determine the byte array size
                    final int expected = readLength - lengthOfSoftClip - (middleOp == 'I' ? lengthOfIndel : (middleOp == 'D' ? -lengthOfIndel : 0));
                    final int indelBasesStart = middleOp != 'M' ? lengthOfFirstM : -1;

                    tests.add(new Object[]{read.getCigar(), expected, middleOp, indelBasesStart, lengthOfIndel});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadToAlignmentByteArrayDataProvider")
    public void testReadToAlignmentByteArrayData(final Cigar cigar, final int expectedLength, final char middleOp, final int startOfIndelBases, final int lengthOfDeletion) {
        final byte[] read = Utils.dupBytes((byte)'A', cigar.getReadLength());
        final byte[] alignment = AlignmentUtils.readToAlignmentByteArray(cigar, read);

        Assert.assertEquals(alignment.length, expectedLength, "Wrong alignment length detected for cigar " + cigar.toString());

        for ( int i = 0; i < alignment.length; i++ ) {
            final byte expectedBase;
            if ( middleOp == 'D' && i >= startOfIndelBases && i < startOfIndelBases + lengthOfDeletion )
                expectedBase = PileupElement.DELETION_BASE;
            else if ( middleOp == 'I' && i == startOfIndelBases - 1 )
                expectedBase = PileupElement.A_FOLLOWED_BY_INSERTION_BASE;
            else
                expectedBase = (byte)'A';
            Assert.assertEquals(alignment[i], expectedBase, "Wrong base detected at position " + i);
        }
    }

    //////////////////////////////////////////
    // Test AlignmentUtils.leftAlignIndel() //
    //////////////////////////////////////////

    @DataProvider(name = "LeftAlignIndelDataProvider")
    public Object[][] makeLeftAlignIndelDataProvider() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final byte[] repeat1Reference = "ABCDEFGHIJKLMNOPXXXXXXXXXXABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat2Reference = "ABCDEFGHIJKLMNOPXYXYXYXYXYABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat3Reference = "ABCDEFGHIJKLMNOPXYZXYZXYZXYZABCDEFGHIJKLMN".getBytes();
        final int referenceLength = repeat1Reference.length;

        for ( int indelStart = 0; indelStart < repeat1Reference.length; indelStart++ ) {
            for ( final int indelSize : Arrays.asList(0, 1, 2, 3, 4) ) {
                for ( final char indelOp : Arrays.asList('D', 'I') ) {

                    if ( indelOp == 'D' && indelStart + indelSize >= repeat1Reference.length )
                        continue;

                    final int readLength = referenceLength - (indelOp == 'D' ? indelSize : -indelSize);

                    // create the original CIGAR string
                    final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
                    read.setCigarString(buildTestCigarString(indelSize == 0 ? 'M' : indelOp, 0, indelStart, indelSize, readLength));
                    final Cigar originalCigar = read.getCigar();

                    final Cigar expectedCigar1 = makeExpectedCigar1(originalCigar, indelOp, indelStart, indelSize, readLength);
                    final byte[] readString1 = makeReadString(repeat1Reference, indelOp, indelStart, indelSize, readLength, 1);
                    tests.add(new Object[]{originalCigar, expectedCigar1, repeat1Reference, readString1, 1});

                    final Cigar expectedCigar2 = makeExpectedCigar2(originalCigar, indelOp, indelStart, indelSize, readLength);
                    final byte[] readString2 = makeReadString(repeat2Reference, indelOp, indelStart, indelSize, readLength, 2);
                    tests.add(new Object[]{originalCigar, expectedCigar2, repeat2Reference, readString2, 2});

                    final Cigar expectedCigar3 = makeExpectedCigar3(originalCigar, indelOp, indelStart, indelSize, readLength);
                    final byte[] readString3 = makeReadString(repeat3Reference, indelOp, indelStart, indelSize, readLength, 3);
                    tests.add(new Object[]{originalCigar, expectedCigar3, repeat3Reference, readString3, 3});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private Cigar makeExpectedCigar1(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength) {
        if ( indelSize == 0 || indelStart < 17 || indelStart > (26 - (indelOp == 'D' ? indelSize : 0)) )
            return originalCigar;

        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
        read.setCigarString(buildTestCigarString(indelOp, 0, 16, indelSize, readLength));
        return read.getCigar();
    }

    private Cigar makeExpectedCigar2(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength) {
        if ( indelStart < 17 || indelStart > (26 - (indelOp == 'D' ? indelSize : 0)) )
            return originalCigar;

        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, readLength);

        if ( indelOp == 'I' && (indelSize == 1 || indelSize == 3) && indelStart % 2 == 1 )
            read.setCigarString(buildTestCigarString(indelOp, 0, Math.max(indelStart - indelSize, 16), indelSize, readLength));
        else if ( (indelSize == 2 || indelSize == 4) && (indelOp == 'D' || indelStart % 2 == 0) )
            read.setCigarString(buildTestCigarString(indelOp, 0, 16, indelSize, readLength));
        else
            return originalCigar;

        return read.getCigar();
    }

    private Cigar makeExpectedCigar3(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength) {
        if ( indelStart < 17 || indelStart > (28 - (indelOp == 'D' ? indelSize : 0)) )
            return originalCigar;

        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, 1, readLength);

        if ( indelSize == 3 && (indelOp == 'D' || indelStart % 3 == 1) )
            read.setCigarString(buildTestCigarString(indelOp, 0, 16, indelSize, readLength));
        else if ( (indelOp == 'I' && indelSize == 4 && indelStart % 3 == 2) ||
                (indelOp == 'I' && indelSize == 2 && indelStart % 3 == 0) ||
                (indelOp == 'I' && indelSize == 1 && indelStart < 28 && indelStart % 3 == 2) )
            read.setCigarString(buildTestCigarString(indelOp, 0, Math.max(indelStart - indelSize, 16), indelSize, readLength));
        else
            return originalCigar;

        return read.getCigar();
    }

    private static byte[] makeReadString(final byte[] reference, final char indelOp, final int indelStart, final int indelSize, final int readLength, final int repeatLength) {
        final byte[] readString = new byte[readLength];

        if ( indelOp == 'D' && indelSize > 0 ) {
            System.arraycopy(reference, 0, readString, 0, indelStart);
            System.arraycopy(reference, indelStart + indelSize, readString, indelStart, readLength - indelStart);
        } else if ( indelOp == 'I' && indelSize > 0 ) {
            System.arraycopy(reference, 0, readString, 0, indelStart);
            for ( int i = 0; i < indelSize; i++ ) {
                if ( i % repeatLength == 0 )
                    readString[indelStart + i] = 'X';
                else if ( i % repeatLength == 1 )
                    readString[indelStart + i] = 'Y';
                else
                    readString[indelStart + i] = 'Z';
            }
            System.arraycopy(reference, indelStart, readString, indelStart + indelSize, readLength - indelStart - indelSize);
        } else {
            System.arraycopy(reference, 0, readString, 0, readLength);
        }

        return readString;
    }

    @Test(dataProvider = "LeftAlignIndelDataProvider", enabled = true)
    public void testLeftAlignIndelData(final Cigar originalCigar, final Cigar expectedCigar, final byte[] reference, final byte[] read, final int repeatLength) {
        final Cigar actualCigar = AlignmentUtils.leftAlignIndel(originalCigar, reference, read, 0, 0, true);
        Assert.assertTrue(expectedCigar.equals(actualCigar), "Wrong left alignment detected for cigar " + originalCigar.toString() + " to " + actualCigar.toString() + " but expected " + expectedCigar.toString() + " with repeat length " + repeatLength);
    }
}
