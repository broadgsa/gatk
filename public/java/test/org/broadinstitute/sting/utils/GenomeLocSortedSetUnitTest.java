package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.Arrays;

/**
 *
 * User: aaron
 * Date: May 22, 2009
 * Time: 2:14:07 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * <p/>
 * Class GenomeLocSetTest
 * <p/>
 * This tests the functions of the GenomeLocSet
 */
public class GenomeLocSortedSetUnitTest extends BaseTest {

    private GenomeLocSortedSet mSortedSet = null;
    private SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    private GenomeLocParser genomeLocParser;
    private String contigOneName;

    @BeforeClass
    public void setup() {
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        contigOneName = header.getSequenceDictionary().getSequence(1).getSequenceName();
    }

    @BeforeMethod
    public void initializeSortedSet() {
        mSortedSet = new GenomeLocSortedSet(genomeLocParser);        
    }

    @Test
    public void testAdd() {
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 0);
        assertTrue(mSortedSet.size() == 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
    }

    @Test
    public void testRemove() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.remove(g);
        assertTrue(mSortedSet.size() == 0);
    }

    @Test
    public void addRegion() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 1, 50);
        mSortedSet.add(g);
        GenomeLoc f = genomeLocParser.createGenomeLoc(contigOneName, 30, 80);
        mSortedSet.addRegion(f);
        assertTrue(mSortedSet.size() == 1);
        
    }


    @Test(expectedExceptions=ReviewedStingException.class)
    public void testAddDuplicate() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.add(g);
    }

    @Test
    public void mergingOverlappingBelow() {
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 0, 50);
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 49, 100);
        assertTrue(mSortedSet.size() == 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.addRegion(e);
        assertTrue(mSortedSet.size() == 1);
        Iterator<GenomeLoc> iter = mSortedSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 0);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == 1);
    }

    @Test
    public void mergingOverlappingAbove() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 0, 50);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 49, 100);
        assertTrue(mSortedSet.size() == 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.addRegion(e);
        assertTrue(mSortedSet.size() == 1);
        Iterator<GenomeLoc> iter = mSortedSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 0);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == 1);
    }

    @Test
    public void deleteAllByRegion() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 1, 100);
        mSortedSet.add(e);
        for (int x = 1; x < 101; x++) {
            GenomeLoc del = genomeLocParser.createGenomeLoc(contigOneName,x,x);
            mSortedSet = mSortedSet.subtractRegions(new GenomeLocSortedSet(genomeLocParser,del));
        }
        assertTrue(mSortedSet.isEmpty());
    }

    @Test
    public void deleteSomeByRegion() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 1, 100);
        mSortedSet.add(e);
        for (int x = 1; x < 50; x++) {
            GenomeLoc del = genomeLocParser.createGenomeLoc(contigOneName,x,x);
            mSortedSet = mSortedSet.subtractRegions(new GenomeLocSortedSet(genomeLocParser,del));
        }
        assertTrue(!mSortedSet.isEmpty());
        assertTrue(mSortedSet.size() == 1);
        GenomeLoc loc = mSortedSet.iterator().next();        
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getStart() == 50);

    }

    @Test
    public void deleteSuperRegion() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 10, 20);
        GenomeLoc g = genomeLocParser.createGenomeLoc(contigOneName, 70, 100);
        mSortedSet.add(g);
        mSortedSet.addRegion(e);
        assertTrue(mSortedSet.size() == 2);
        // now delete a region
        GenomeLoc d = genomeLocParser.createGenomeLoc(contigOneName, 15, 75);
        mSortedSet = mSortedSet.subtractRegions(new GenomeLocSortedSet(genomeLocParser,d));
        Iterator<GenomeLoc> iter = mSortedSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 10);
        assertTrue(loc.getStop() == 14);
        assertTrue(loc.getContigIndex() == 1);

        loc = iter.next();
        assertTrue(loc.getStart() == 76);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == 1);
    }

    @Test
    public void substractComplexExample() {
        GenomeLoc e = genomeLocParser.createGenomeLoc(contigOneName, 1, 20);
        mSortedSet.add(e);

        GenomeLoc r1 = genomeLocParser.createGenomeLoc(contigOneName, 3, 5);
        GenomeLoc r2 = genomeLocParser.createGenomeLoc(contigOneName, 10, 12);
        GenomeLoc r3 = genomeLocParser.createGenomeLoc(contigOneName, 16, 18);
        GenomeLocSortedSet toExclude = new GenomeLocSortedSet(genomeLocParser,Arrays.asList(r1, r2, r3));

        GenomeLocSortedSet remaining = mSortedSet.subtractRegions(toExclude);
//        logger.debug("Initial   " + mSortedSet);
//        logger.debug("Exclude   " + toExclude);
//        logger.debug("Remaining " + remaining);

        assertEquals(mSortedSet.coveredSize(), 20);
        assertEquals(toExclude.coveredSize(), 9);
        assertEquals(remaining.coveredSize(), 11);

        Iterator<GenomeLoc> it = remaining.iterator();
        GenomeLoc p1 = it.next();
        GenomeLoc p2 = it.next();
        GenomeLoc p3 = it.next();
        GenomeLoc p4 = it.next();

        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 1, 2), p1);
        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 6, 9), p2);
        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 13, 15), p3);
        assertEquals(genomeLocParser.createGenomeLoc(contigOneName, 19, 20), p4);
    }

    private void testSizeBeforeLocX(int pos, int size) {
        GenomeLoc test = genomeLocParser.createGenomeLoc(contigOneName, pos, pos);
        assertEquals(mSortedSet.sizeBeforeLoc(test), size, String.format("X pos=%d size=%d", pos, size));
    }

    @Test
    public void testSizeBeforeLoc() {
        GenomeLoc r1 = genomeLocParser.createGenomeLoc(contigOneName, 3, 5);
        GenomeLoc r2 = genomeLocParser.createGenomeLoc(contigOneName, 10, 12);
        GenomeLoc r3 = genomeLocParser.createGenomeLoc(contigOneName, 16, 18);
        mSortedSet.addAll(Arrays.asList(r1,r2,r3));

        testSizeBeforeLocX(2, 0);
        testSizeBeforeLocX(3, 0);
        testSizeBeforeLocX(4, 1);
        testSizeBeforeLocX(5, 2);
        testSizeBeforeLocX(6, 3);

        testSizeBeforeLocX(10, 3);
        testSizeBeforeLocX(11, 4);
        testSizeBeforeLocX(12, 5);
        testSizeBeforeLocX(13, 6);
        testSizeBeforeLocX(15, 6);

        testSizeBeforeLocX(16, 6);
        testSizeBeforeLocX(17, 7);
        testSizeBeforeLocX(18, 8);
        testSizeBeforeLocX(19, 9);
        testSizeBeforeLocX(50, 9);
        testSizeBeforeLocX(50, (int)mSortedSet.coveredSize());
    }


    @Test
    public void fromSequenceDictionary() {
        mSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(this.header.getSequenceDictionary());
        // we should have sequence
        assertTrue(mSortedSet.size() == GenomeLocSortedSetUnitTest.NUMBER_OF_CHROMOSOMES);
        int seqNumber = 0;
        for (GenomeLoc loc : mSortedSet) {
            assertTrue(loc.getStart() == 1);
            assertTrue(loc.getStop() == GenomeLocSortedSetUnitTest.CHROMOSOME_SIZE);
            assertTrue(loc.getContigIndex() == seqNumber);
            ++seqNumber;
        }
        assertTrue(seqNumber == GenomeLocSortedSetUnitTest.NUMBER_OF_CHROMOSOMES);
    }
}
