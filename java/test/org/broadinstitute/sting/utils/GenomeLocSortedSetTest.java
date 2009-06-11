package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.Test;

import java.util.Iterator;

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
public class GenomeLocSortedSetTest extends BaseTest {

    private GenomeLocSortedSet mSortedSet = null;
    private SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    @Before
    public void setup() {
        GenomeLoc.setupRefContigOrdering(header.getSequenceDictionary());
        mSortedSet = new GenomeLocSortedSet();
    }

    @Test
    public void testAdd() {
        GenomeLoc g = new GenomeLoc(1, 0, 0);
        assertTrue(mSortedSet.size() == 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
    }

    @Test
    public void testRemove() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = new GenomeLoc(1, 0, 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.remove(g);
        assertTrue(mSortedSet.size() == 0);
    }

    @Test
    public void addRegion() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = new GenomeLoc(1, 1, 50);
        mSortedSet.add(g);
        GenomeLoc f = new GenomeLoc(1, 30, 80);
        mSortedSet.addRegion(f);
        assertTrue(mSortedSet.size() == 1);
        
    }


    @Test(expected = IllegalArgumentException.class)
    public void testAddDupplicate() {
        assertTrue(mSortedSet.size() == 0);
        GenomeLoc g = new GenomeLoc(1, 0, 0);
        mSortedSet.add(g);
        assertTrue(mSortedSet.size() == 1);
        mSortedSet.add(g);
    }

    @Test
    public void mergingOverlappingBelow() {
        GenomeLoc g = new GenomeLoc(1, 0, 50);
        GenomeLoc e = new GenomeLoc(1, 49, 100);
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
        GenomeLoc e = new GenomeLoc(1, 0, 50);
        GenomeLoc g = new GenomeLoc(1, 49, 100);
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
    public void deleteSubRegion() {
        GenomeLoc e = new GenomeLoc(1, 0, 50);
        GenomeLoc g = new GenomeLoc(1, 49, 100);
        mSortedSet.add(g);
        mSortedSet.addRegion(e);

        // now delete a region
        GenomeLoc d = new GenomeLoc(1, 25, 75);
        mSortedSet.removeRegion(d);
        Iterator<GenomeLoc> iter = mSortedSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 0);
        assertTrue(loc.getStop() == 24);
        assertTrue(loc.getContigIndex() == 1);

        loc = iter.next();
        assertTrue(loc.getStart() == 76);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == 1);
    }

    @Test
    public void deleteAllByRegion() {
        GenomeLoc e = new GenomeLoc(1, 1, 100);
        mSortedSet.add(e);
        for (int x = 1; x < 101; x++) {
            GenomeLoc del = new GenomeLoc(1,x,x);
            mSortedSet.removeRegion(del);
        }
        assertTrue(mSortedSet.isEmpty());
    }
    @Test
    public void deleteSomeByRegion() {
        GenomeLoc e = new GenomeLoc(1, 1, 100);
        mSortedSet.add(e);
        for (int x = 1; x < 50; x++) {
            GenomeLoc del = new GenomeLoc(1,x,x);
            mSortedSet.removeRegion(del);
        }
        assertTrue(!mSortedSet.isEmpty());
        assertTrue(mSortedSet.size() == 1);
        GenomeLoc loc = mSortedSet.iterator().next();        
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getStart() == 50);

    }

    @Test
    public void deleteSuperRegion() {
        GenomeLoc e = new GenomeLoc(1, 10, 20);
        GenomeLoc g = new GenomeLoc(1, 70, 100);
        mSortedSet.add(g);
        mSortedSet.addRegion(e);
        assertTrue(mSortedSet.size() == 2);
        // now delete a region
        GenomeLoc d = new GenomeLoc(1, 15, 75);
        mSortedSet.removeRegion(d);
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
    public void fromSequenceDictionary() {
        mSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(this.header.getSequenceDictionary());
        // we should have sequence
        assertTrue(mSortedSet.size() == GenomeLocSortedSetTest.NUMBER_OF_CHROMOSOMES);
        int seqNumber = 0;
        for (GenomeLoc loc : mSortedSet) {
            assertTrue(loc.getStart() == 1);
            assertTrue(loc.getStop() == GenomeLocSortedSetTest.CHROMOSOME_SIZE);
            assertTrue(loc.getContigIndex() == seqNumber);
            ++seqNumber;
        }
        assertTrue(seqNumber == GenomeLocSortedSetTest.NUMBER_OF_CHROMOSOMES);
    }
}
