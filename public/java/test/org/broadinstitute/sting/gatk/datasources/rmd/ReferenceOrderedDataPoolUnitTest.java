package org.broadinstitute.sting.gatk.datasources.rmd;

import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet.RMDStorageType;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;

import static org.testng.Assert.assertTrue;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.picard.reference.IndexedFastaSequenceFile;
/**
 * User: hanna
 * Date: May 21, 2009
 * Time: 11:03:04 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Test the contents and number of iterators in the pool.
 */

public class ReferenceOrderedDataPoolUnitTest extends BaseTest {

    private RMDTriplet triplet = null;
    private RMDTrackBuilder builder = null;

    private IndexedFastaSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    private GenomeLoc testSite1;
    private GenomeLoc testSite2;
    private GenomeLoc testSite3;

    private GenomeLoc testInterval1; // an interval matching testSite1 -> testSite2 for queries
    private GenomeLoc testInterval2; // an interval matching testSite2 -> testSite3 for queries


    @BeforeClass
    public void init() throws FileNotFoundException {
        seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq);

        testSite1 = genomeLocParser.createGenomeLoc("chrM",10);
        testSite2 = genomeLocParser.createGenomeLoc("chrM",20);
        testSite3 = genomeLocParser.createGenomeLoc("chrM",30);
        testInterval1 = genomeLocParser.createGenomeLoc("chrM",10,20);
        testInterval2 = genomeLocParser.createGenomeLoc("chrM",20,30);
    }

    @BeforeMethod
    public void setUp() {
        String fileName = testDir + "TabularDataTest.dat";

        // check to see if we have an index, if so, delete it
        File indexFileName = new File(testDir + "TabularDataTest.dat.idx");
        if (indexFileName.exists()) indexFileName.delete();  

        triplet = new RMDTriplet("tableTest","Table",fileName,RMDStorageType.FILE,new Tags());
        builder = new RMDTrackBuilder(seq.getSequenceDictionary(),genomeLocParser,null);
    }

    @Test
    public void testCreateSingleIterator() {
        ResourcePool iteratorPool = new ReferenceOrderedDataPool(triplet,builder,seq.getSequenceDictionary(),genomeLocParser,false);
        LocationAwareSeekableRODIterator iterator = (LocationAwareSeekableRODIterator)iteratorPool.iterator( new MappedStreamSegment(testSite1) );

        Assert.assertEquals(iteratorPool.numIterators(), 1, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 0, "Number of available iterators in the pool is incorrect");

        TableFeature datum = (TableFeature)iterator.next().get(0).getUnderlyingObject();

        assertTrue(datum.getLocation().equals(testSite1));
        assertTrue(datum.get("COL1").equals("A"));
        assertTrue(datum.get("COL2").equals("B"));
        assertTrue(datum.get("COL3").equals("C"));

        iteratorPool.release(iterator);

        Assert.assertEquals(iteratorPool.numIterators(), 1, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 1, "Number of available iterators in the pool is incorrect");
    }

    @Test
    public void testCreateMultipleIterators() {
        ReferenceOrderedQueryDataPool iteratorPool = new ReferenceOrderedQueryDataPool(triplet,builder,seq.getSequenceDictionary(),genomeLocParser);
        LocationAwareSeekableRODIterator iterator1 = iteratorPool.iterator( new MappedStreamSegment(testInterval1) );

        // Create a new iterator at position 2.
        LocationAwareSeekableRODIterator iterator2 = iteratorPool.iterator( new MappedStreamSegment(testInterval2) );

        Assert.assertEquals(iteratorPool.numIterators(), 2, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 0, "Number of available iterators in the pool is incorrect");

        // Test out-of-order access: first iterator2, then iterator1.
        // Ugh...first call to a region needs to be a seek.
        TableFeature datum = (TableFeature)iterator2.seekForward(testSite2).get(0).getUnderlyingObject();
        assertTrue(datum.getLocation().equals(testSite2));
        assertTrue(datum.get("COL1").equals("C"));
        assertTrue(datum.get("COL2").equals("D"));
        assertTrue(datum.get("COL3").equals("E"));

        datum = (TableFeature)iterator1.next().get(0).getUnderlyingObject();
        assertTrue(datum.getLocation().equals(testSite1));
        assertTrue(datum.get("COL1").equals("A"));
        assertTrue(datum.get("COL2").equals("B"));
        assertTrue(datum.get("COL3").equals("C"));

        // Advance iterator2, and make sure both iterator's contents are still correct.
        datum = (TableFeature)iterator2.next().get(0).getUnderlyingObject();
        assertTrue(datum.getLocation().equals(testSite3));
        assertTrue(datum.get("COL1").equals("F"));
        assertTrue(datum.get("COL2").equals("G"));
        assertTrue(datum.get("COL3").equals("H"));

        datum = (TableFeature)iterator1.next().get(0).getUnderlyingObject();
        assertTrue(datum.getLocation().equals(testSite2));
        assertTrue(datum.get("COL1").equals("C"));
        assertTrue(datum.get("COL2").equals("D"));
        assertTrue(datum.get("COL3").equals("E"));

        // Cleanup, and make sure the number of iterators dies appropriately.
        iteratorPool.release(iterator1);

        Assert.assertEquals(iteratorPool.numIterators(), 2, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 1, "Number of available iterators in the pool is incorrect");

        iteratorPool.release(iterator2);

        Assert.assertEquals(iteratorPool.numIterators(), 2, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 2, "Number of available iterators in the pool is incorrect");
    }

    @Test
    public void testIteratorConservation() {
        ReferenceOrderedDataPool iteratorPool = new ReferenceOrderedDataPool(triplet,builder,seq.getSequenceDictionary(),genomeLocParser,false);
        LocationAwareSeekableRODIterator iterator = iteratorPool.iterator( new MappedStreamSegment(testSite1) );

        Assert.assertEquals(iteratorPool.numIterators(), 1, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 0, "Number of available iterators in the pool is incorrect");

        TableFeature datum = (TableFeature)iterator.next().get(0).getUnderlyingObject();
        assertTrue(datum.getLocation().equals(testSite1));
        assertTrue(datum.get("COL1").equals("A"));
        assertTrue(datum.get("COL2").equals("B"));
        assertTrue(datum.get("COL3").equals("C"));

        iteratorPool.release(iterator);

        // Create another iterator after the current iterator.
        iterator = iteratorPool.iterator( new MappedStreamSegment(testSite3) );

        // Make sure that the previously acquired iterator was reused.
        Assert.assertEquals(iteratorPool.numIterators(), 1, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 0, "Number of available iterators in the pool is incorrect");

        datum = (TableFeature)iterator.seekForward(testSite3).get(0).getUnderlyingObject();
        assertTrue(datum.getLocation().equals(testSite3));
        assertTrue(datum.get("COL1").equals("F"));
        assertTrue(datum.get("COL2").equals("G"));
        assertTrue(datum.get("COL3").equals("H"));

        iteratorPool.release(iterator);

        Assert.assertEquals(iteratorPool.numIterators(), 1, "Number of iterators in the pool is incorrect");
        Assert.assertEquals(iteratorPool.numAvailableIterators(), 1, "Number of available iterators in the pool is incorrect");
    }
}
