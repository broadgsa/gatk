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

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.GATKChunk;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;

/**
 *
 */
public class FilePointerUnitTest extends BaseTest {
    private IndexedFastaSequenceFile seq;
    private GenomeLocParser genomeLocParser;
    private SAMReaderID readerID = new SAMReaderID("samFile",new Tags());

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @BeforeMethod
    public void doForEachTest() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq.getSequenceDictionary());
    }

    @Test
    public void testFilePointerCombineDisjoint() {
        FilePointer one = new FilePointer(genomeLocParser.createGenomeLoc("chr1",1,5));
        one.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(0,1)));
        FilePointer two = new FilePointer(genomeLocParser.createGenomeLoc("chr1",6,10));
        two.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(1,2)));

        FilePointer result = new FilePointer(genomeLocParser.createGenomeLoc("chr1",1,10));
        result.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(0,2)));

        Assert.assertEquals(one.combine(genomeLocParser,two),result,"Combination of two file pointers is incorrect");
        Assert.assertEquals(two.combine(genomeLocParser,one),result,"Combination of two file pointers is incorrect");
    }

    @Test
    public void testFilePointerCombineJoint() {
        FilePointer one = new FilePointer(genomeLocParser.createGenomeLoc("chr1",1,5));
        one.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(0,2)));
        FilePointer two = new FilePointer(genomeLocParser.createGenomeLoc("chr1",2,6));
        two.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(1,3)));

        FilePointer result = new FilePointer(genomeLocParser.createGenomeLoc("chr1",1,6));
        result.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(0,3)));        

        Assert.assertEquals(one.combine(genomeLocParser,two),result,"Combination of two file pointers is incorrect");
        Assert.assertEquals(two.combine(genomeLocParser,one),result,"Combination of two file pointers is incorrect");
    }

    @Test
    public void testFilePointerCombineOneSided() {
        FilePointer filePointer = new FilePointer(genomeLocParser.createGenomeLoc("chr1",1,5));
        filePointer.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(0,1)));
        FilePointer empty = new FilePointer(genomeLocParser.createGenomeLoc("chr1",6,10));
        // Do not add file spans to empty result

        FilePointer result = new FilePointer(genomeLocParser.createGenomeLoc("chr1",1,10));
        result.addFileSpans(readerID,new GATKBAMFileSpan(new GATKChunk(0,1)));
        Assert.assertEquals(filePointer.combine(genomeLocParser,empty),result,"Combination of two file pointers is incorrect");
        Assert.assertEquals(empty.combine(genomeLocParser,filePointer),result,"Combination of two file pointers is incorrect");
    }
}
