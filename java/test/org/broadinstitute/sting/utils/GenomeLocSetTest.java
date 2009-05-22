package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSamUtils;
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
 * @date May 22, 2009
 * <p/>
 * Class GenomeLocSetTest
 * <p/>
 * This tests the functions of the GenomeLocSet
 */
public class GenomeLocSetTest extends BaseTest {

    private GenomeLocSet mSet = null;
    private SAMFileHeader header = ArtificialSamUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    @Before
    public void setup() {
        GenomeLoc.setupRefContigOrdering(header.getSequenceDictionary());
        mSet = new GenomeLocSet();
    }

    @Test
    public void testAdd() {
        GenomeLoc g = new GenomeLoc(STARTING_CHROMOSOME, 0, 0);
        assertTrue(mSet.size() == 0);
        mSet.add(g);
        assertTrue(mSet.size() == STARTING_CHROMOSOME);
    }

    @Test
    public void testRemove() {
        assertTrue(mSet.size() == 0);
        GenomeLoc g = new GenomeLoc(STARTING_CHROMOSOME, 0, 0);
        mSet.add(g);
        assertTrue(mSet.size() == STARTING_CHROMOSOME);
        mSet.remove(g);
        assertTrue(mSet.size() == 0);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testAddDupplicate() {
        assertTrue(mSet.size() == 0);
        GenomeLoc g = new GenomeLoc(STARTING_CHROMOSOME, 0, 0);
        mSet.add(g);
        assertTrue(mSet.size() == STARTING_CHROMOSOME);
        mSet.add(g);
    }

    @Test
    public void mergingOverlappingBelow() {
        GenomeLoc g = new GenomeLoc(STARTING_CHROMOSOME, 0, 50);
        GenomeLoc e = new GenomeLoc(STARTING_CHROMOSOME, 49, 100);
        assertTrue(mSet.size() == 0);
        mSet.add(g);
        assertTrue(mSet.size() == STARTING_CHROMOSOME);
        mSet.addRegion(e);
        assertTrue(mSet.size() == STARTING_CHROMOSOME);
        Iterator<GenomeLoc> iter = mSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 0);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == STARTING_CHROMOSOME);
    }

    @Test
    public void mergingOverlappingAbove() {
        GenomeLoc e = new GenomeLoc(STARTING_CHROMOSOME, 0, 50);
        GenomeLoc g = new GenomeLoc(STARTING_CHROMOSOME, 49, 100);
        assertTrue(mSet.size() == 0);
        mSet.add(g);
        assertTrue(mSet.size() == STARTING_CHROMOSOME);
        mSet.addRegion(e);
        assertTrue(mSet.size() == STARTING_CHROMOSOME);
        Iterator<GenomeLoc> iter = mSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 0);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == STARTING_CHROMOSOME);
    }

    @Test
    public void deleteSubRegion() {
        GenomeLoc e = new GenomeLoc(STARTING_CHROMOSOME, 0, 50);
        GenomeLoc g = new GenomeLoc(STARTING_CHROMOSOME, 49, 100);
        mSet.add(g);
        mSet.addRegion(e);

        // now delete a region
        GenomeLoc d = new GenomeLoc(STARTING_CHROMOSOME, 25, 75);
        mSet.removeRegion(d);
        Iterator<GenomeLoc> iter = mSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 0);
        assertTrue(loc.getStop() == 24);
        assertTrue(loc.getContigIndex() == STARTING_CHROMOSOME);

        loc = iter.next();
        assertTrue(loc.getStart() == 76);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == STARTING_CHROMOSOME);
    }
    @Test
    public void deleteSuperRegion() {
        GenomeLoc e = new GenomeLoc(STARTING_CHROMOSOME, 10, 20);
        GenomeLoc g = new GenomeLoc(STARTING_CHROMOSOME, 70, 100);
        mSet.add(g);
        mSet.addRegion(e);
        assertTrue(mSet.size() == 2);
        // now delete a region
        GenomeLoc d = new GenomeLoc(STARTING_CHROMOSOME, 15, 75);
        mSet.removeRegion(d);
        Iterator<GenomeLoc> iter = mSet.iterator();
        GenomeLoc loc = iter.next();
        assertTrue(loc.getStart() == 10);
        assertTrue(loc.getStop() == 14);
        assertTrue(loc.getContigIndex() == STARTING_CHROMOSOME);

        loc = iter.next();
        assertTrue(loc.getStart() == 76);
        assertTrue(loc.getStop() == 100);
        assertTrue(loc.getContigIndex() == STARTING_CHROMOSOME);
    }

    @Test
    public void fromSequenceDictionary() {
        mSet = GenomeLocSet.createSetFromSequenceDictionary(this.header.getSequenceDictionary());
        // we should have sequence
        assertTrue(mSet.size() == GenomeLocSetTest.NUMBER_OF_CHROMOSOMES);
        int seqNumber = 0;
        for (GenomeLoc loc : mSet) {
            assertTrue(loc.getStart() == 1);
            assertTrue(loc.getStop() == GenomeLocSetTest.CHROMOSOME_SIZE);
            assertTrue(loc.getContigIndex() == seqNumber);
            ++seqNumber;
        }
        assertTrue(seqNumber == GenomeLocSetTest.NUMBER_OF_CHROMOSOMES);
    }
}
