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

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IntervalOverlapFilteringIteratorUnitTest {

    private SAMFileHeader header;
    private GenomeLoc firstContig;
    private GenomeLoc secondContig;

    /** Basic aligned and mapped read. */
    private SAMRecord readMapped;

    /** Read with no contig specified in the read, -L UNMAPPED */
    private SAMRecord readNoReference;

    /** This read has a start position, but is flagged that it's not mapped. */
    private SAMRecord readUnmappedFlag;

    /** This read is from the second contig. */
    private SAMRecord readSecondContig;

    /** This read says it's aligned, but actually has an unknown start. */
    private SAMRecord readUnknownStart;

    /** The above reads in the order one would expect to find them in a sorted BAM. */
    private List<SAMRecord> testReads;

    @BeforeClass
    public void init() {
        header = ArtificialSAMUtils.createArtificialSamHeader(3, 1, ArtificialSAMUtils.DEFAULT_READ_LENGTH * 2);
        GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        SAMSequenceRecord record;

        record = header.getSequence(0);
        firstContig = genomeLocParser.createGenomeLoc(record.getSequenceName(), 1, record.getSequenceLength());
        record = header.getSequence(1);
        secondContig = genomeLocParser.createGenomeLoc(record.getSequenceName(), 1, record.getSequenceLength());

        readMapped = createMappedRead("mapped", 1);

        readUnmappedFlag = createMappedRead("unmappedFlagged", 2);
        readUnmappedFlag.setReadUnmappedFlag(true);

        readSecondContig = createMappedRead("secondContig", 3);
        readSecondContig.setReferenceName(secondContig.getContig());

        /* This read says it's aligned, but to a contig not in the header. */
        SAMRecord readUnknownContig = createMappedRead("unknownContig", 4);
        readUnknownContig.setReferenceName("unknownContig");

        readUnknownStart = createMappedRead("unknownStart", 1);
        readUnknownStart.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);

        readNoReference = createUnmappedRead("unmappedNoReference");

        testReads = new ArrayList<SAMRecord>();
        testReads.add(readMapped);
        testReads.add(readUnmappedFlag);
        testReads.add(readUnknownStart);
        testReads.add(readSecondContig);
        testReads.add(readUnknownContig);
        testReads.add(readNoReference);
    }

    @DataProvider(name = "filteringIteratorTestData")
    public Object[][] getFilteringIteratorTestData() {
        return new Object[][] {
                new Object[] {Arrays.asList(firstContig), Arrays.asList(readMapped, readUnmappedFlag, readUnknownStart)},
                new Object[] {Arrays.asList(GenomeLoc.UNMAPPED), Arrays.asList(readNoReference)},
                new Object[] {Arrays.asList(firstContig, secondContig), Arrays.asList(readMapped, readUnmappedFlag, readUnknownStart, readSecondContig)}
        };
    }

    @Test(dataProvider = "filteringIteratorTestData")
    public void testFilteringIterator(List<GenomeLoc> locs, List<SAMRecord> expected) {
        IntervalOverlapFilteringIterator filterIter = new IntervalOverlapFilteringIterator(
                ArtificialSAMUtils.createReadIterator(testReads), locs);

        List<SAMRecord> actual = new ArrayList<SAMRecord>();
        while (filterIter.hasNext()) {
            actual.add(filterIter.next());
        }
        Assert.assertEquals(actual, expected);
    }

    @Test(expectedExceptions = ReviewedGATKException.class)
    public void testMappedAndUnmapped() {
        new IntervalOverlapFilteringIterator(
                ArtificialSAMUtils.createReadIterator(testReads),
                Arrays.asList(firstContig, GenomeLoc.UNMAPPED));
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
