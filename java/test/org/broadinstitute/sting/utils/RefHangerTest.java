// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Basic unit test for RefHanger
 *
 * interface functions:
    public RefHanger();
    public void clear();
    protected int getLeftOffset();
    protected int getRightOffset();
    protected int getOffset(GenomeLoc loc);
    public GenomeLoc getLeftLoc();
    public GenomeLoc getRightLoc();
    public boolean hasLocation(GenomeLoc loc);
    public boolean isEmpty();
    public boolean hasHangers();
    public Hanger popLeft();
    public void dropLeft();
    public Hanger getLeft();
    public Hanger getHanger(int relativePos);
    public Hanger getHanger(GenomeLoc pos);
    public int size();
    public void pushLeft(GenomeLoc pos);
    public void pushLeft(GenomeLoc pos, T datum);
    public void pushLeft(GenomeLoc pos, ArrayList<T> data);
    public void pushRight(GenomeLoc pos);
    public void pushRight(GenomeLoc pos, T datum);
    public void pushRight(GenomeLoc pos, ArrayList<T> data);
    public boolean ensurePos(GenomeLoc pos);
    public void extendData(GenomeLoc pos, T datum);
    public void addData(List<GenomeLoc> positions, List<T> dataByPos);
    public void expandingPut1(final GenomeLoc loc, T datum);
    public void printState();
    public void expandingPut(GenomeLoc pos, T datum);

 */
public class RefHangerTest extends BaseTest {
    private static FastaSequenceFile2 seq;
    private GenomeLoc startLoc;
    private RefHanger<Integer> emptyHanger;

    /**
     * chrM:1 is the start.
     * 1: 1 2 3
     * 2: 4 5
     * 3: 6
     * 4: 7 8
     * 5: 9 10
     */
    private static RefHanger<Integer> filledHanger;
    private static List<Integer> l1, l2, l3, l4, l5;
    private static GenomeLoc p1, p2, p3, p4, p5;

    @BeforeClass
    public static void init() {
        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq);

        System.out.printf("Filled hanger is %n%s%n", makeFilledHanger());
    }

    private static RefHanger makeFilledHanger() {
        RefHanger<Integer> filledHanger = new RefHanger<Integer>();
        l1 = Arrays.asList(1, 2, 3);
        l2 = Arrays.asList(4, 5);
        l3 = Arrays.asList(6);
        l4 = Arrays.asList(7, 8);
        l5 = Arrays.asList(9, 10);
        p1 = GenomeLocParser.createGenomeLoc(0, 1, 1);
        p2 = new GenomeLoc(p1).nextLoc();
        p3 = new GenomeLoc(p2).nextLoc();
        p4 = new GenomeLoc(p3).nextLoc();
        p5 = new GenomeLoc(p4).nextLoc();

        filledHanger.addDataList(Arrays.asList(p1, p2, p3, p4, p5),
                                 Arrays.asList(l1, l2, l3, l4, l5));
        return filledHanger;
    }

    @Before
    public void setupHanger() {
        startLoc = GenomeLocParser.createGenomeLoc(0, 1, 1);  // chrM 1
        emptyHanger = new RefHanger<Integer>();
        filledHanger = makeFilledHanger();

        //System.out.printf("Filled hanger is %n%s%n", filledHanger);

        //GenomeLoc two = new GenomeLoc();
    }

    /**
     * Tests that we got a string parameter in correctly
     */
    @Test
    public void testEmpty() {
        logger.warn("Executing testEmpty");
        assertTrue(emptyHanger.size() == 0);
        assertTrue(emptyHanger.getLeftOffset() == 0);
        assertTrue(emptyHanger.getRightOffset() == 0);
    }

    @Test (expected=IndexOutOfBoundsException.class )
    public void testEmptyGet() {
        logger.warn("Executing testEmptyGet");
        emptyHanger.getHanger(0);
    }

    @Test
    public void testFilledFilling() {
        logger.warn("Executing testFilledFilling");
        testBaseFilledHanger();
    }

    private void testBaseFilledHanger() {
        assertTrue(filledHanger.size() == 5);
        assertTrue(! filledHanger.isEmpty());
        assertTrue(filledHanger.hasHangers());
        assertTrue(filledHanger.getLeftOffset() == 0);
        assertTrue(filledHanger.getRightOffset() == 4);
        assertTrue(filledHanger.getOffset(p1) == 0);
        assertTrue(filledHanger.getOffset(p2) == 1);
        assertTrue(filledHanger.getOffset(p3) == 2);
        assertTrue(filledHanger.getOffset(p4) == 3);
        assertTrue(filledHanger.getOffset(p5) == 4);
        assertTrue(filledHanger.getLeftLoc().compareTo(p1) == 0);
        assertTrue(filledHanger.getRightLoc().compareTo(p5) == 0);
        assertTrue(filledHanger.getLoc(2).compareTo(p3) == 0);
        assertTrue(filledHanger.hasLocation(p1));
        assertTrue(filledHanger.hasLocation(p2));
        assertTrue(filledHanger.hasLocation(p3));
        assertTrue(filledHanger.hasLocation(p4));
        assertTrue(filledHanger.hasLocation(p5));
        assertTrue(! filledHanger.hasLocation(GenomeLocParser.createGenomeLoc(0, 6, 6)));

        assertTrue(filledHanger.getHanger(0) != null);
        assertTrue(filledHanger.getHanger(1) != null);
        assertTrue((int)(Integer)filledHanger.getHanger(0).data.get(0) == 1);
        assertTrue((int)(Integer)filledHanger.getHanger(1).data.get(0) == 4);
        assertTrue((int)(Integer)filledHanger.getHanger(p1).data.get(0) == 1);
        assertTrue((int)(Integer)filledHanger.getHanger(p2).data.get(0) == 4);
    }

    @Test
    public void testClear() {
        logger.warn("Executing testClear");        
        assertTrue(filledHanger.size() == 5);
        assertTrue(! filledHanger.isEmpty());
        assertTrue(filledHanger.hasHangers());

        filledHanger.clear();

        assertTrue(filledHanger.size() == 0);
        assertTrue(filledHanger.isEmpty());
        assertTrue(! filledHanger.hasHangers());
    }

    @Test
    public void testPopLeft() {
        logger.warn("Executing testPopLeft");
        filledHanger.popLeft();
        assertTrue(filledHanger.size() == 4);
        assertTrue(filledHanger.getRightOffset() == 3);
        assertTrue(filledHanger.getOffset(p2) == 0);
        assertTrue(filledHanger.getOffset(p3) == 1);
        assertTrue(filledHanger.getOffset(p4) == 2);
        assertTrue(filledHanger.getOffset(p5) == 3);
        assertTrue(! filledHanger.hasLocation(p1));
        assertTrue(filledHanger.hasLocation(p2));
        assertTrue(filledHanger.hasLocation(p3));
        assertTrue(filledHanger.hasLocation(p4));
        assertTrue(filledHanger.hasLocation(p5));
    }

    @Test
    public void testDoublePopLeft() {
        logger.warn("Executing testDoublePopLeft");
        filledHanger.popLeft();
        filledHanger.popLeft();
        assertTrue(filledHanger.size() == 3);
        assertTrue(filledHanger.getRightOffset() == 2);
        assertTrue(filledHanger.getOffset(p3) == 0);
        assertTrue(filledHanger.getOffset(p4) == 1);
        assertTrue(filledHanger.getOffset(p5) == 2);
        assertTrue(! filledHanger.hasLocation(p1));
        assertTrue(! filledHanger.hasLocation(p2));
        assertTrue(filledHanger.hasLocation(p3));
        assertTrue(filledHanger.hasLocation(p4));
        assertTrue(filledHanger.hasLocation(p5));
    }

    @Test
    public void testPopPushLeft() {
        logger.warn("Executing testPopPushLeft");
        filledHanger.popLeft();
        filledHanger.pushLeft(p1, l1);
        testBaseFilledHanger();
    }
}