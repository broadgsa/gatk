/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.utils.pileup;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Test routines for read-backed pileup.
 */
public class ReadBackedPileupTest {
    /**
     * Ensure that basic read group splitting works.
     */
    @Test
    public void testSplitByReadGroup() {
        SAMReadGroupRecord readGroupOne = new SAMReadGroupRecord("rg1");
        SAMReadGroupRecord readGroupTwo = new SAMReadGroupRecord("rg2");

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
        header.addReadGroup(readGroupOne);
        header.addReadGroup(readGroupTwo);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);
        read1.setAttribute("RG",readGroupOne.getId());
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,10);
        read2.setAttribute("RG",readGroupTwo.getId());
        SAMRecord read3 = ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,10);
        read3.setAttribute("RG",readGroupOne.getId());
        SAMRecord read4 = ArtificialSAMUtils.createArtificialRead(header,"read4",0,1,10);
        read4.setAttribute("RG",readGroupTwo.getId());
        SAMRecord read5 = ArtificialSAMUtils.createArtificialRead(header,"read5",0,1,10);
        read5.setAttribute("RG",readGroupTwo.getId());
        SAMRecord read6 = ArtificialSAMUtils.createArtificialRead(header,"read6",0,1,10);
        read6.setAttribute("RG",readGroupOne.getId());
        SAMRecord read7 = ArtificialSAMUtils.createArtificialRead(header,"read7",0,1,10);
        read7.setAttribute("RG",readGroupOne.getId());

        ReadBackedPileup pileup = new ReadBackedPileupImpl(null,
                                                           Arrays.asList(read1,read2,read3,read4,read5,read6,read7),
                                                           Arrays.asList(1,1,1,1,1,1,1));

        ReadBackedPileup rg1Pileup = pileup.getPileupForReadGroup("rg1");
        List<SAMRecord> rg1Reads = rg1Pileup.getReads();
        Assert.assertEquals("Wrong number of reads in read group rg1",4,rg1Reads.size());
        Assert.assertEquals("Read " + read1.getReadName() + " should be in rg1 but isn't",read1,rg1Reads.get(0));
        Assert.assertEquals("Read " + read3.getReadName() + " should be in rg1 but isn't",read3,rg1Reads.get(1));
        Assert.assertEquals("Read " + read6.getReadName() + " should be in rg1 but isn't",read6,rg1Reads.get(2));
        Assert.assertEquals("Read " + read7.getReadName() + " should be in rg1 but isn't",read7,rg1Reads.get(3));

        ReadBackedPileup rg2Pileup = pileup.getPileupForReadGroup("rg2");
        List<SAMRecord> rg2Reads = rg2Pileup.getReads();        
        Assert.assertEquals("Wrong number of reads in read group rg2",3,rg2Reads.size());
        Assert.assertEquals("Read " + read2.getReadName() + " should be in rg2 but isn't",read2,rg2Reads.get(0));
        Assert.assertEquals("Read " + read4.getReadName() + " should be in rg2 but isn't",read4,rg2Reads.get(1));
        Assert.assertEquals("Read " + read5.getReadName() + " should be in rg2 but isn't",read5,rg2Reads.get(2));
    }

    /**
     * Ensure that splitting read groups still works when dealing with null read groups.
     */
    @Test
    public void testSplitByNullReadGroups() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,10);
        SAMRecord read3 = ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,10);

        ReadBackedPileup pileup = new ReadBackedPileupImpl(null,
                                                           Arrays.asList(read1,read2,read3),
                                                           Arrays.asList(1,1,1));

        ReadBackedPileup nullRgPileup = pileup.getPileupForReadGroup(null);
        List<SAMRecord> nullRgReads = nullRgPileup.getReads();
        Assert.assertEquals("Wrong number of reads in null read group",3,nullRgPileup.size());
        Assert.assertEquals("Read " + read1.getReadName() + " should be in null rg but isn't",read1,nullRgReads.get(0));
        Assert.assertEquals("Read " + read2.getReadName() + " should be in null rg but isn't",read2,nullRgReads.get(1));
        Assert.assertEquals("Read " + read3.getReadName() + " should be in null rg but isn't",read3,nullRgReads.get(2));

        ReadBackedPileup rg1Pileup = pileup.getPileupForReadGroup("rg1");
        Assert.assertNull("Pileup for non-existent read group should return null",rg1Pileup);
    }

    /**
     * Ensure that splitting read groups still works when dealing with a sample-split pileup.
     */
    @Test
    public void testSplitBySample() {
        SAMReadGroupRecord readGroupOne = new SAMReadGroupRecord("rg1");
        readGroupOne.setSample("sample1");
        SAMReadGroupRecord readGroupTwo = new SAMReadGroupRecord("rg2");
        readGroupTwo.setSample("sample2");

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
        header.addReadGroup(readGroupOne);
        header.addReadGroup(readGroupTwo);

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,10);
        read1.setAttribute("RG",readGroupOne.getId());
        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,10);
        read2.setAttribute("RG",readGroupTwo.getId());
        SAMRecord read3 = ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,10);
        read3.setAttribute("RG",readGroupOne.getId());
        SAMRecord read4 = ArtificialSAMUtils.createArtificialRead(header,"read4",0,1,10);
        read4.setAttribute("RG",readGroupTwo.getId());

        ReadBackedPileupImpl sample1Pileup = new ReadBackedPileupImpl(null,
                                                                      Arrays.asList(read1,read3),
                                                                      Arrays.asList(1,1));
        ReadBackedPileupImpl sample2Pileup = new ReadBackedPileupImpl(null,
                                                                      Arrays.asList(read2,read4),
                                                                      Arrays.asList(1,1));
        Map<String,ReadBackedPileupImpl> sampleToPileupMap = new HashMap<String,ReadBackedPileupImpl>();
        sampleToPileupMap.put(readGroupOne.getSample(),sample1Pileup);
        sampleToPileupMap.put(readGroupTwo.getSample(),sample2Pileup);

        ReadBackedPileup compositePileup = new ReadBackedPileupImpl(null,sampleToPileupMap);

        ReadBackedPileup rg1Pileup = compositePileup.getPileupForReadGroup("rg1");
        List<SAMRecord> rg1Reads = rg1Pileup.getReads();

        Assert.assertEquals("Wrong number of reads in read group rg1",2,rg1Reads.size());
        Assert.assertEquals("Read " + read1.getReadName() + " should be in rg1 but isn't",read1,rg1Reads.get(0));
        Assert.assertEquals("Read " + read3.getReadName() + " should be in rg1 but isn't",read3,rg1Reads.get(1));

        ReadBackedPileup rg2Pileup = compositePileup.getPileupForReadGroup("rg2");
        List<SAMRecord> rg2Reads = rg2Pileup.getReads();

        Assert.assertEquals("Wrong number of reads in read group rg2",2,rg1Reads.size());
        Assert.assertEquals("Read " + read2.getReadName() + " should be in rg2 but isn't",read2,rg2Reads.get(0));
        Assert.assertEquals("Read " + read4.getReadName() + " should be in rg2 but isn't",read4,rg2Reads.get(1));
    }
}
