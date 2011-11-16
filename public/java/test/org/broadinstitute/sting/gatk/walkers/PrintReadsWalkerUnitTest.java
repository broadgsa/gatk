package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.sam.ArtificialReadsTraversal;
import org.broadinstitute.sting.utils.sam.ArtificialSAMFileWriter;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class PrintReadsWalkerUnitTest
 *         <p/>
 *         This tests the print reads walker, using the artificial reads traversal
 */
public class PrintReadsWalkerUnitTest extends BaseTest {

    /**
     * our private fake reads traversal.  This traversal seeds the
     * the walker with a specified number of fake reads, which are named sequentially
     */
    private ArtificialReadsTraversal trav;
    private int readTotal = 0;
    //private char bases[] = {'a', 't'};
    private ReferenceContext bases = null;
    //private ReferenceContext ref = new ReferenceContext()

    @BeforeMethod
    public void before() {
        trav = new ArtificialReadsTraversal();
        readTotal = ( ( trav.endingChr - trav.startingChr ) + 1 ) * trav.readsPerChr + trav.unMappedReads;
    }

    /** test that we get out the same number of reads we put in */
    @Test
    public void testReadCount() {
        PrintReadsWalker walker = new PrintReadsWalker();
        ArtificialSAMFileWriter writer = new ArtificialSAMFileWriter();
        walker.out = writer;

        trav.traverse(walker, null, writer);
        assertEquals(writer.getRecords().size(), readTotal);
    }

    /** test that we're ok with a null read */
    @Test
    public void testNullRead() {
        PrintReadsWalker walker = new PrintReadsWalker();
        ArtificialSAMFileWriter writer = new ArtificialSAMFileWriter();
        walker.out = writer;

        SAMRecord rec = walker.map(bases, null, null);
        assertTrue(rec == null);
    }

    /** tes that we get the read we put into the map function */
    @Test
    public void testReturnRead() {
        PrintReadsWalker walker = new PrintReadsWalker();
        ArtificialSAMFileWriter writer = new ArtificialSAMFileWriter();
        walker.out = writer;

        SAMFileHeader head = ArtificialSAMUtils.createArtificialSamHeader(3,1,1000);
        GATKSAMRecord rec = ArtificialSAMUtils.createArtificialRead(head, "FakeRead", 1, 1, 50);
        SAMRecord ret = walker.map(bases, rec, null);
        assertTrue(ret == rec);
        assertTrue(ret.getReadName().equals(rec.getReadName()));
    }

    /** test that the read makes it to the output source */
    @Test
    public void testReducingRead() {
        PrintReadsWalker walker = new PrintReadsWalker();
        ArtificialSAMFileWriter writer = new ArtificialSAMFileWriter();
        walker.out = writer;

        SAMFileHeader head = ArtificialSAMUtils.createArtificialSamHeader(3,1,1000);
        SAMRecord rec = ArtificialSAMUtils.createArtificialRead(head, "FakeRead", 1, 1, 50);
        SAMRecord ret = walker.map(bases, null,null);
        walker.reduce(ret,writer);

        assertTrue(writer.getRecords().size() == 1);
    }

    
}
