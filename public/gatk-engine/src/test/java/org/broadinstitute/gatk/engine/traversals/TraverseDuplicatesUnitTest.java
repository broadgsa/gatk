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

package org.broadinstitute.gatk.engine.traversals;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.testng.Assert;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;

import org.testng.annotations.BeforeMethod;

import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;


/**
 * @author aaron
 *         <p/>
 *         Class TraverseDuplicatesUnitTest
 *         <p/>
 *         test the meat of the traverse dupplicates.
 */
public class TraverseDuplicatesUnitTest extends BaseTest {

    private TraverseDuplicates obj = new TraverseDuplicates();
    private SAMFileHeader header;
    private GenomeLocParser genomeLocParser;
    private GenomeAnalysisEngine engine;
    private File refFile = new File(validationDataLocation + "Homo_sapiens_assembly17.fasta");


    @BeforeMethod
    public void doBefore() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser =new GenomeLocParser(header.getSequenceDictionary());

        engine = new GenomeAnalysisEngine();
        engine.setReferenceDataSource(refFile);
        engine.setGenomeLocParser(genomeLocParser);
        
        obj.initialize(engine, null);
    }

    @Test
    public void testAllDuplicatesNoPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ" + x, 0, 1, 100);
            read.setDuplicateReadFlag(true);
            list.add(read);
        }
        Set<List<SAMRecord>> myPairings = obj.uniqueReadSets(list);
        Assert.assertEquals(myPairings.size(), 1);
        Assert.assertEquals(myPairings.iterator().next().size(), 10); // dup's
    }

    @Test
    public void testNoDuplicatesNoPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ" + x, 0, 1, 100);
            read.setDuplicateReadFlag(false);
            list.add(read);
        }

        Set<List<SAMRecord>> myPairing = obj.uniqueReadSets(list);
        Assert.assertEquals(myPairing.size(), 10); // unique
    }

    @Test
    public void testFiftyFiftyNoPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 5; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ" + x, 0, 1, 100);
            read.setDuplicateReadFlag(true);
            list.add(read);
        }
        for (int x = 10; x < 15; x++)
            list.add(ArtificialSAMUtils.createArtificialRead(header, String.valueOf(x), 0, x, 100));

        Set<List<SAMRecord>> myPairing = obj.uniqueReadSets(list);
        Assert.assertEquals(myPairing.size(), 6);  // unique
    }

    @Test
    public void testAllDuplicatesAllPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ"+ x, 0, 1, 100);
            read.setDuplicateReadFlag(true);
            read.setMateAlignmentStart(100);
            read.setMateReferenceIndex(0);
            read.setReadPairedFlag(true);
            list.add(read);
        }

        Set<List<SAMRecord>> myPairing = obj.uniqueReadSets(list);
        Assert.assertEquals(myPairing.size(), 1);  // unique
    }

    @Test
    public void testNoDuplicatesAllPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ"+ x, 0, 1, 100);
            if (x == 0) read.setDuplicateReadFlag(true); // one is a dup but (next line)
            read.setMateAlignmentStart(100); // they all have a shared start and mate start so they're dup's
            read.setMateReferenceIndex(0);
            read.setReadPairedFlag(true);
            list.add(read);
        }

        Set<List<SAMRecord>> myPairing = obj.uniqueReadSets(list);
        Assert.assertEquals(myPairing.size(), 1);  // unique
    }

    @Test
    public void testAllDuplicatesAllPairsDifferentPairedEnd() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ" + x, 0, 1, 100);
            if (x == 0) read.setDuplicateReadFlag(true); // one is a dup
            read.setMateAlignmentStart(100 + x);
            read.setMateReferenceIndex(0);
            read.setReadPairedFlag(true);
            list.add(read);
        }

        Set<List<SAMRecord>> myPairing = obj.uniqueReadSets(list);
        Assert.assertEquals(myPairing.size(), 10);  // unique
    }
}
