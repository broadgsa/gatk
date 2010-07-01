// our package
package org.broadinstitute.sting.gatk.refdata;


// the imports for unit testing.

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeatureIterator;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import static org.junit.Assert.assertTrue;

/**
 * Basic unit test for TabularROD
 *
 */
public class TabularRODUnitTest extends BaseTest {
    private static ReferenceSequenceFile seq;
    private ReferenceOrderedData ROD;
    private LocationAwareSeekableRODIterator iter;


    @BeforeClass
    public static void init() throws FileNotFoundException {
        // sequence
        seq = new IndexedFastaSequenceFile(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    @Before
    public void setupTabularROD() {
        TabularROD.setDelimiter(TabularROD.DEFAULT_DELIMITER, TabularROD.DEFAULT_DELIMITER_REGEX);        
        File file = new File(testDir + "TabularDataTest.dat");
        ROD = new ReferenceOrderedData("tableTest", file, TabularROD.class);
        iter = new SeekableRODIterator(new GATKFeatureIterator(ROD.iterator()));

    }

    @Test
    public void test1() {
        logger.warn("Executing test1");
        RODRecordList oneList = iter.next();
        TabularROD one = (TabularROD)oneList.get(0).getUnderlyingObject();
        assertTrue(one.size() == 4);
        assertTrue(one.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 10)));
        assertTrue(one.get("COL1").equals("A"));
        assertTrue(one.get("COL2").equals("B"));        
        assertTrue(one.get("COL3").equals("C"));
    }

    @Test
    public void test2() {
        logger.warn("Executing test2");
        RODRecordList oneList = iter.next();
        RODRecordList twoList = iter.next();
        TabularROD one = (TabularROD)oneList.get(0).getUnderlyingObject();
        TabularROD two = (TabularROD)twoList.get(0).getUnderlyingObject();
        assertTrue(two.size() == 4);
        assertTrue(two.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 20)));
        assertTrue(two.get("COL1").equals("C"));
        assertTrue(two.get("COL2").equals("D"));
        assertTrue(two.get("COL3").equals("E"));
    }

    @Test
    public void test3() {
        logger.warn("Executing test3");
        RODRecordList oneList = iter.next();
        RODRecordList twoList = iter.next();
        RODRecordList threeList = iter.next();
        TabularROD one = (TabularROD)oneList.get(0).getUnderlyingObject();
        TabularROD two = (TabularROD)twoList.get(0).getUnderlyingObject();
        TabularROD three = (TabularROD)threeList.get(0).getUnderlyingObject();
        assertTrue(three.size() == 4);
        assertTrue(three.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 30)));
        assertTrue(three.get("COL1").equals("F"));
        assertTrue(three.get("COL2").equals("G"));
        assertTrue(three.get("COL3").equals("H"));
    }

    @Test
    public void testDone() {
        logger.warn("Executing testDone");
        RODRecordList oneList = iter.next();
        RODRecordList twoList = iter.next();
        RODRecordList threeList = iter.next();
        TabularROD one = (TabularROD)oneList.get(0).getUnderlyingObject();
        TabularROD two = (TabularROD)twoList.get(0).getUnderlyingObject();
        TabularROD three = (TabularROD)threeList.get(0).getUnderlyingObject();
        assertTrue(!iter.hasNext());
    }

    @Test
    public void testSeek() {
        logger.warn("Executing testSeek");
        RODRecordList twoList = iter.seekForward(GenomeLocParser.createGenomeLoc("chrM", 20));
        TabularROD two = (TabularROD)twoList.get(0).getUnderlyingObject();
        assertTrue(two.size() == 4);
        assertTrue(two.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 20)));
        assertTrue(two.get("COL1").equals("C"));
        assertTrue(two.get("COL2").equals("D"));
        assertTrue(two.get("COL3").equals("E"));
    }

    @Test
    public void testToString() {
        logger.warn("Executing testToString");
        RODRecordList oneList = iter.next();
        TabularROD one = (TabularROD)oneList.get(0).getUnderlyingObject();
        assertTrue(one.toString().equals("chrM:10\tA\tB\tC"));
    }

    // Didn't change the delimiter
    @Test (expected = RuntimeException.class)
    public void testDelim1() {
        File file2 = new File(testDir + "TabularDataTest2.dat");        
        ReferenceOrderedData ROD_commas = new ReferenceOrderedData("tableTest", file2, TabularROD.class);
        LocationAwareSeekableRODIterator iter_commas = new SeekableRODIterator(new GATKFeatureIterator(ROD_commas.iterator()));

        logger.warn("Executing testDelim1");
        RODRecordList one2List = iter_commas.next();
        TabularROD one2 = (TabularROD)one2List.get(0).getUnderlyingObject();
        assertTrue(one2.size() == 5);
        assertTrue(one2.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 10)));
        assertTrue(one2.get("COL1").equals("A"));
        assertTrue(one2.get("COL2").equals("B"));
        assertTrue(one2.get("COL3").equals("C"));
        assertTrue(one2.get("COL4").equals("1"));
    }

    @Test
    public void testDelim2() {
        TabularROD.setDelimiter(",",",");
        File file2 = new File(testDir + "TabularDataTest2.dat");
        ReferenceOrderedData ROD_commas = new ReferenceOrderedData("tableTest", file2, TabularROD.class);
        LocationAwareSeekableRODIterator iter_commas = new SeekableRODIterator(new GATKFeatureIterator(ROD_commas.iterator()));

        logger.warn("Executing testDelim1");
        RODRecordList one2List = iter_commas.next();
        TabularROD one2 = (TabularROD)one2List.get(0).getUnderlyingObject();
        assertTrue(one2.size() == 5);
        assertTrue(one2.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 10)));
        assertTrue(one2.get("COL1").equals("A"));
        assertTrue(one2.get("COL2").equals("B"));
        assertTrue(one2.get("COL3").equals("C"));
        assertTrue(one2.get("COL4").equals("1"));
    }

    @Test
    public void testCreation() {
        logger.warn("Executing testCreation");
        ArrayList<String> header = new ArrayList<String>(Arrays.asList("HEADER", "col1", "col2", "col3"));
        assertTrue(TabularROD.headerString(header).equals("HEADER\tcol1\tcol2\tcol3"));
        String rowData = String.format("%d %d %d", 1, 2, 3);
        TabularROD row = new TabularROD("myName", header, GenomeLocParser.createGenomeLoc("chrM", 1), rowData.split(" "));
        System.out.println(">>>>> " + row.toString());
        assertTrue(row.toString().equals("chrM:1\t1\t2\t3"));
    }

    @Test
    public void testCreationAndWriting() throws FileNotFoundException {
        logger.warn("Executing testCreationAndWriting");

        File outputFile = new File(testDir + "testTabularRodOutputTemp.dat");
        PrintStream out = new PrintStream(new FileOutputStream(outputFile));

        ArrayList<String> header = new ArrayList<String>(Arrays.asList("HEADER", "col1", "col2", "col3"));
        out.println(TabularROD.commentString("Hello, created from test"));
        out.println(TabularROD.commentString(""));
        out.println(TabularROD.headerString(header));

        String rowData = String.format("%d %d %d", 1, 2, 3);
        TabularROD row = new TabularROD("myName", header, GenomeLocParser.createGenomeLoc("chrM", 1), rowData.split(" "));
        out.println(row.toString());

        rowData = String.format("%d %d %d", 3, 4, 5);
        row = new TabularROD("myName", header, GenomeLocParser.createGenomeLoc("chrM", 2), rowData.split(" "));
        out.println(row.toString());

        ReferenceOrderedData ROD_commas = new ReferenceOrderedData("tableTest", outputFile, TabularROD.class);
        LocationAwareSeekableRODIterator iter_commas = new SeekableRODIterator(new GATKFeatureIterator(ROD_commas.iterator()));

        RODRecordList oneList = iter_commas.next();
        TabularROD one = (TabularROD)oneList.get(0).getUnderlyingObject();
        assertTrue(one.size() == 4);
        assertTrue(one.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 1)));
        assertTrue(one.get("col1").equals("1"));
        assertTrue(one.get("col2").equals("2"));
        assertTrue(one.get("col3").equals("3"));

        RODRecordList twoList = iter_commas.next();
        TabularROD two = (TabularROD)twoList.get(0).getUnderlyingObject();
        assertTrue(two.size() == 4);
        assertTrue(two.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 2)));
        assertTrue(two.get("col1").equals("3"));
        assertTrue(two.get("col2").equals("4"));
        assertTrue(two.get("col3").equals("5"));
    }

/*    @Test (expected=RuntimeException.class )
    public void testBadHeader1() {
        logger.warn("Executing testBadHeader1");
        ArrayList<String> header = new ArrayList<String>();
        TabularROD row = new TabularROD("myName", header, GenomeLocParser.createGenomeLoc("chrM", 1));
    }*/

/*    @Test (expected=RuntimeException.class )
    public void testBadHeader2() {
        logger.warn("Executing testBadHeader2");
        ArrayList<String> header = new ArrayList<String>(Arrays.asList("col1", "col2", "col3"));
        TabularROD row = new TabularROD("myName", header, GenomeLocParser.createGenomeLoc("chrM", 1));
    }*/

    @Test (expected=RuntimeException.class )
    public void testBadData1() {
        logger.warn("Executing testBadData1");
        ArrayList<String> header = new ArrayList<String>(Arrays.asList("HEADER", "col1", "col2", "col3"));
        assertTrue(TabularROD.headerString(header).equals("HEADER\tcol1\tcol2\tcol3"));
        String rowData = String.format("%d %d %d %d", 1, 2, 3, 4);
        TabularROD row = new TabularROD("myName", header, GenomeLocParser.createGenomeLoc("chrM", 1), rowData.split(" "));
    }
}