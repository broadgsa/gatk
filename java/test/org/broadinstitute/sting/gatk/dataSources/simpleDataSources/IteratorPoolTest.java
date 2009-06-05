package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Assert;
import static org.junit.Assert.assertTrue;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.TabularROD;
import org.broadinstitute.sting.gatk.refdata.RODIterator;

import java.io.File;
import java.io.FileNotFoundException;
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

public class IteratorPoolTest extends BaseTest {

    private static File sequenceFile = new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta");

    private ReferenceOrderedData<? extends ReferenceOrderedDatum> rod = null;

    private final GenomeLoc testSite1 = new GenomeLoc("chrM",10);
    private final GenomeLoc testSite2 = new GenomeLoc("chrM",20);
    private final GenomeLoc testSite3 = new GenomeLoc("chrM",30);

    @BeforeClass
    public static void init() throws FileNotFoundException {
        GenomeLoc.setupRefContigOrdering(new IndexedFastaSequenceFile(sequenceFile));
        TabularROD.setDelimiter(TabularROD.DEFAULT_DELIMITER, TabularROD.DEFAULT_DELIMITER_REGEX);
    }

    @Before
    public void setUp() {
        File file = new File(testDir + "TabularDataTest.dat");
        rod = new ReferenceOrderedData("tableTest", file, TabularROD.class);
    }

    @Test
    public void testCreateSingleIterator() {
        IteratorPool iteratorPool = new IteratorPool(rod); 
        RODIterator iterator = (RODIterator)iteratorPool.iterator( testSite1 );

        Assert.assertEquals("Number of iterators in the pool is incorrect", 1, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 0, iteratorPool.numAvailableIterators());

        TabularROD datum = (TabularROD)iterator.next();

        assertTrue(datum.getLocation().equals(testSite1));
        assertTrue(datum.get("COL1").equals("A"));
        assertTrue(datum.get("COL2").equals("B"));
        assertTrue(datum.get("COL3").equals("C"));

        iteratorPool.close(iterator);

        Assert.assertEquals("Number of iterators in the pool is incorrect", 1, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 1, iteratorPool.numAvailableIterators());        
    }

    @Test
    public void testCreateMultipleIterators() {
        IteratorPool iteratorPool = new IteratorPool(rod);
        RODIterator iterator1 = (RODIterator)iteratorPool.iterator( testSite1 );

        // Create a new iterator at position 2.
        RODIterator iterator2 = iteratorPool.iterator( testSite2 );

        Assert.assertEquals("Number of iterators in the pool is incorrect", 2, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 0, iteratorPool.numAvailableIterators());

        // Test out-of-order access: first iterator2, then iterator1.
        // Ugh...first call to a region needs to be a seek. 
        TabularROD datum = (TabularROD)iterator2.seekForward(testSite2);
        assertTrue(datum.getLocation().equals(testSite2));
        assertTrue(datum.get("COL1").equals("C"));
        assertTrue(datum.get("COL2").equals("D"));
        assertTrue(datum.get("COL3").equals("E"));

        datum = (TabularROD)iterator1.next();
        assertTrue(datum.getLocation().equals(testSite1));
        assertTrue(datum.get("COL1").equals("A"));
        assertTrue(datum.get("COL2").equals("B"));
        assertTrue(datum.get("COL3").equals("C"));

        // Advance iterator2, and make sure both iterator's contents are still correct.
        datum = (TabularROD)iterator2.next();
        assertTrue(datum.getLocation().equals(testSite3));
        assertTrue(datum.get("COL1").equals("F"));
        assertTrue(datum.get("COL2").equals("G"));
        assertTrue(datum.get("COL3").equals("H"));

        datum = (TabularROD)iterator1.next();
        assertTrue(datum.getLocation().equals(testSite2));
        assertTrue(datum.get("COL1").equals("C"));
        assertTrue(datum.get("COL2").equals("D"));
        assertTrue(datum.get("COL3").equals("E"));

        // Cleanup, and make sure the number of iterators dies appropriately.
        iteratorPool.close(iterator1);

        Assert.assertEquals("Number of iterators in the pool is incorrect", 2, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 1, iteratorPool.numAvailableIterators());

        iteratorPool.close(iterator2);

        Assert.assertEquals("Number of iterators in the pool is incorrect", 2, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 2, iteratorPool.numAvailableIterators());        
    }

    @Test
    public void testIteratorConservation() {
        IteratorPool iteratorPool = new IteratorPool(rod);
        RODIterator iterator = (RODIterator)iteratorPool.iterator( testSite1 );

        Assert.assertEquals("Number of iterators in the pool is incorrect", 1, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 0, iteratorPool.numAvailableIterators());

        TabularROD datum = (TabularROD)iterator.next();
        assertTrue(datum.getLocation().equals(testSite1));
        assertTrue(datum.get("COL1").equals("A"));
        assertTrue(datum.get("COL2").equals("B"));
        assertTrue(datum.get("COL3").equals("C"));

        iteratorPool.close(iterator);

        // Create another iterator after the current iterator.
        iterator = iteratorPool.iterator(testSite3);

        // Make sure that the previously acquired iterator was reused.
        Assert.assertEquals("Number of iterators in the pool is incorrect", 1, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 0, iteratorPool.numAvailableIterators());

        datum = (TabularROD)iterator.seekForward(testSite3);
        assertTrue(datum.getLocation().equals(testSite3));
        assertTrue(datum.get("COL1").equals("F"));
        assertTrue(datum.get("COL2").equals("G"));
        assertTrue(datum.get("COL3").equals("H"));

        iteratorPool.close(iterator);

        Assert.assertEquals("Number of iterators in the pool is incorrect", 1, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 1, iteratorPool.numAvailableIterators());
    }

    @Test
    public void testIteratorCreation() {
        IteratorPool iteratorPool = new IteratorPool(rod);
        RODIterator iterator = (RODIterator)iteratorPool.iterator( testSite3 );

        Assert.assertEquals("Number of iterators in the pool is incorrect", 1, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 0, iteratorPool.numAvailableIterators());

        TabularROD datum = (TabularROD)iterator.seekForward(testSite3);
        assertTrue(datum.getLocation().equals(testSite3));
        assertTrue(datum.get("COL1").equals("F"));
        assertTrue(datum.get("COL2").equals("G"));
        assertTrue(datum.get("COL3").equals("H"));

        iteratorPool.close(iterator);

        // Create another iterator after the current iterator.
        iterator = iteratorPool.iterator(testSite1);

        // Make sure that the previously acquired iterator was reused.
        Assert.assertEquals("Number of iterators in the pool is incorrect", 2, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 1, iteratorPool.numAvailableIterators());

        datum = (TabularROD)iterator.next();
        assertTrue(datum.getLocation().equals(testSite1));
        assertTrue(datum.get("COL1").equals("A"));
        assertTrue(datum.get("COL2").equals("B"));
        assertTrue(datum.get("COL3").equals("C"));

        iteratorPool.close(iterator);

        Assert.assertEquals("Number of iterators in the pool is incorrect", 2, iteratorPool.numIterators());
        Assert.assertEquals("Number of available iterators in the pool is incorrect", 2, iteratorPool.numAvailableIterators());
    }

}
