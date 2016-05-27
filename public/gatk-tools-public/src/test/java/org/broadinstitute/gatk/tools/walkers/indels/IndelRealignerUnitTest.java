/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;

public class IndelRealignerUnitTest extends BaseTest {

    private SAMFileHeader header;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        final IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
    }

    @Test
    public void realignAtContigBorderTest() {
        final int contigEnd = header.getSequence(0).getSequenceLength();
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "goodRead", 0, contigEnd - 1, 2);
        read.setCigarString("2M");
        Assert.assertEquals(IndelRealigner.realignmentProducesBadAlignment(read, contigEnd), false);
        read.setCigarString("1M1D1M");
        Assert.assertEquals(IndelRealigner.realignmentProducesBadAlignment(read, contigEnd), true);
    }

}
