/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Mar 2, 2011
 * Time: 9:48:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class VerifyingSamIteratorUnitTest {
    private SAMFileHeader samFileHeader;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() {
        SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary();
        sequenceDictionary.addSequence(new SAMSequenceRecord("1",500));
        sequenceDictionary.addSequence(new SAMSequenceRecord("2",500));

        samFileHeader = new SAMFileHeader();
        samFileHeader.setSequenceDictionary(sequenceDictionary);

        genomeLocParser = new GenomeLocParser(sequenceDictionary);
    }

    @Test
    public void testSortedReadsBasic() {
        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read1",getContig(0).getSequenceIndex(),1,10);
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read2",getContig(0).getSequenceIndex(),2,10);
        List<SAMRecord> reads = Arrays.asList(read1,read2);

        VerifyingSamIterator iterator = new VerifyingSamIterator(genomeLocParser,StingSAMIteratorAdapter.adapt(reads.iterator()));

        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");
        Assert.assertSame(iterator.next(),read1,"Incorrect read in read 1 position");
        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");
        Assert.assertSame(iterator.next(),read2,"Incorrect read in read 2 position");
        Assert.assertFalse(iterator.hasNext(),"Too many reads in iterator");
    }

    @Test
    public void testSortedReadsAcrossContigs() {
        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read1",getContig(0).getSequenceIndex(),2,10);
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read2",getContig(1).getSequenceIndex(),1,10);
        List<SAMRecord> reads = Arrays.asList(read1,read2);

        VerifyingSamIterator iterator = new VerifyingSamIterator(genomeLocParser,StingSAMIteratorAdapter.adapt(reads.iterator()));

        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");
        Assert.assertSame(iterator.next(),read1,"Incorrect read in read 1 position");
        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");
        Assert.assertSame(iterator.next(),read2,"Incorrect read in read 2 position");
        Assert.assertFalse(iterator.hasNext(),"Too many reads in iterator");
    }

    @Test(expectedExceptions=UserException.MissortedBAM.class)
    public void testImproperlySortedReads() {
        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read1",getContig(0).getSequenceIndex(),2,10);
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read2",getContig(0).getSequenceIndex(),1,10);
        List<SAMRecord> reads = Arrays.asList(read1,read2);

        VerifyingSamIterator iterator = new VerifyingSamIterator(genomeLocParser,StingSAMIteratorAdapter.adapt(reads.iterator()));

        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");
        Assert.assertSame(iterator.next(),read1,"Incorrect read in read 1 position");
        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");

        // Should trigger MissortedBAM exception.
        iterator.next();
    }

    @Test(expectedExceptions=UserException.MalformedBAM.class)
    public void testInvalidAlignment() {
        // Create an invalid alignment state.
        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read1",getContig(0).getSequenceIndex(),1,10);
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(samFileHeader,"read1",getContig(0).getSequenceIndex(),2,10);
        read1.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        List<SAMRecord> reads = Arrays.asList(read1,read2);

        VerifyingSamIterator iterator = new VerifyingSamIterator(genomeLocParser,StingSAMIteratorAdapter.adapt(reads.iterator()));

        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");
        Assert.assertSame(iterator.next(),read1,"Incorrect read in read 1 position");
        Assert.assertTrue(iterator.hasNext(),"Insufficient reads");

        // Should trigger MalformedBAM exception.
        iterator.next();
    }

    private SAMSequenceRecord getContig(final int contigIndex) {
        return samFileHeader.getSequence(contigIndex);            
    }
}
