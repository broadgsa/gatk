// our package
package org.broadinstitute.sting.gatk.refdata;


// the imports for unit testing.

import org.junit.*;
import static org.junit.Assert.assertTrue;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.utils.RefHanger;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Basic unit test for TabularROD
 *
 */
public class TabularRODTest extends BaseTest {
    private static FastaSequenceFile2 seq;
    private ReferenceOrderedData ROD;
    private ReferenceOrderedData.RODIterator iter;
    private ReferenceOrderedData ROD2;
    private ReferenceOrderedData.RODIterator iter2;

    @BeforeClass
    public static void init() {
        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLoc.setupRefContigOrdering(seq);
    }

    @Before
    public void setupTabularROD() {
        File file = new File(testDir + "TabularDataTest.dat");
        ROD = new ReferenceOrderedData("tableTest", file, TabularROD.class);
        iter = ROD.iterator();
        
        File file2 = new File(testDir + "TabularDataTest2.dat");
        ROD2 = new ReferenceOrderedData("tableTest", file2, TabularROD.class);
        iter2 = ROD2.iterator();
    }

    @Test
    public void test1() {
        logger.warn("Executing test1");
        TabularROD one = (TabularROD)iter.next();
        assertTrue(one.size() == 3);
        assertTrue(one.getLocation().equals(new GenomeLoc("chrM", 10)));
        assertTrue(one.get("COL1").equals("A"));
        assertTrue(one.get("COL2").equals("B"));        
        assertTrue(one.get("COL3").equals("C"));
    }

    @Test
    public void test2() {
        logger.warn("Executing test2");
        TabularROD one = (TabularROD)iter.next();
        TabularROD two = (TabularROD)iter.next();
        assertTrue(two.size() == 3);
        assertTrue(two.getLocation().equals(new GenomeLoc("chrM", 20)));
        assertTrue(two.get("COL1").equals("C"));
        assertTrue(two.get("COL2").equals("D"));
        assertTrue(two.get("COL3").equals("E"));
    }

    @Test
    public void test3() {
        logger.warn("Executing test3");
        TabularROD one = (TabularROD)iter.next();
        TabularROD two = (TabularROD)iter.next();
        TabularROD three = (TabularROD)iter.next();
        assertTrue(three.size() == 3);
        assertTrue(three.getLocation().equals(new GenomeLoc("chrM", 30)));
        assertTrue(three.get("COL1").equals("F"));
        assertTrue(three.get("COL2").equals("G"));
        assertTrue(three.get("COL3").equals("H"));
    }

    @Test
    public void testDone() {
        logger.warn("Executing testDone");
        TabularROD one = (TabularROD)iter.next();
        TabularROD two = (TabularROD)iter.next();
        TabularROD three = (TabularROD)iter.next();
        assertTrue(!iter.hasNext());
    }

    @Test
    public void testSeek() {
        logger.warn("Executing testSeek");
        TabularROD two = (TabularROD)iter.seekForward(new GenomeLoc("chrM", 20));
        assertTrue(two.size() == 3);
        assertTrue(two.getLocation().equals(new GenomeLoc("chrM", 20)));
        assertTrue(two.get("COL1").equals("C"));
        assertTrue(two.get("COL2").equals("D"));
        assertTrue(two.get("COL3").equals("E"));
    }

    @Test
    public void testToString() {
        logger.warn("Executing testToString");
        TabularROD one = (TabularROD)iter.next();
        assertTrue(one.toString().equals("chrM:10\tA\tB\tC"));
    }

    @Test
    public void test2p1() {
        logger.warn("Executing test2p1");
        TabularROD one2 = (TabularROD)iter2.next();
        assertTrue(one2.size() == 4);
        assertTrue(one2.getLocation().equals(new GenomeLoc("chrM", 10)));
        assertTrue(one2.get("COL1").equals("A"));
        assertTrue(one2.get("COL2").equals("B"));
        assertTrue(one2.get("COL3").equals("C"));
        assertTrue(one2.get("COL4").equals("1"));
    }
}