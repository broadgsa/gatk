/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.sam;

import junit.framework.Assert;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

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
}
