// our package
package org.broadinstitute.sting.gatk.refdata;


// the imports for unit testing.

import org.junit.*;
import static org.junit.Assert.assertTrue;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.io.File;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.ArrayList;

/**
 * Basic unit test for TabularROD
 *
 */
public class TabularRODTest extends BaseTest {
    private static FastaSequenceFile2 seq;
    private ReferenceOrderedData ROD;
    private RODIterator iter;


    @BeforeClass
    public static void init() {
        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    @Before
    public void setupTabularROD() {
        TabularROD.setDelimiter(TabularROD.DEFAULT_DELIMITER, TabularROD.DEFAULT_DELIMITER_REGEX);        
        File file = new File(testDir + "TabularDataTest.dat");
        ROD = new ReferenceOrderedData("tableTest", file, TabularROD.class);
        iter = ROD.iterator();

    }

    @Test
    public void test1() {
        logger.warn("Executing test1");
        TabularROD one = (TabularROD)iter.next();
        assertTrue(one.size() == 4);
        assertTrue(one.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 10)));
        assertTrue(one.get("COL1").equals("A"));
        assertTrue(one.get("COL2").equals("B"));        
        assertTrue(one.get("COL3").equals("C"));
    }

    @Test
    public void test2() {
        logger.warn("Executing test2");
        TabularROD one = (TabularROD)iter.next();
        TabularROD two = (TabularROD)iter.next();
        assertTrue(two.size() == 4);
        assertTrue(two.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 20)));
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
        assertTrue(three.size() == 4);
        assertTrue(three.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 30)));
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
        TabularROD two = (TabularROD)iter.seekForward(GenomeLocParser.createGenomeLoc("chrM", 20));
        assertTrue(two.size() == 4);
        assertTrue(two.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 20)));
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

    // Didn't change the delimiter
    @Test (expected = RuntimeException.class)
    public void testDelim1() {
        File file2 = new File(testDir + "TabularDataTest2.dat");        
        ReferenceOrderedData ROD_commas = new ReferenceOrderedData("tableTest", file2, TabularROD.class);
        RODIterator iter_commas = ROD_commas.iterator();

        logger.warn("Executing testDelim1");
        TabularROD one2 = (TabularROD)iter_commas.next();
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
        RODIterator iter_commas = ROD_commas.iterator();

        logger.warn("Executing testDelim1");
        TabularROD one2 = (TabularROD)iter_commas.next();
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
        RODIterator iter_commas = ROD_commas.iterator();

        TabularROD one = (TabularROD)iter_commas.next();
        assertTrue(one.size() == 4);
        assertTrue(one.getLocation().equals(GenomeLocParser.createGenomeLoc("chrM", 1)));
        assertTrue(one.get("col1").equals("1"));
        assertTrue(one.get("col2").equals("2"));
        assertTrue(one.get("col3").equals("3"));

        TabularROD two = (TabularROD)iter_commas.next();
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