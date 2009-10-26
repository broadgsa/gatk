package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class TraverseDuplicatesTest
 *         <p/>
 *         test the meat of the traverse dupplicates.
 */
public class TraverseDuplicatesTest extends BaseTest {

    private TraverseDuplicates obj = new TraverseDuplicates();
    private SAMFileHeader header;


    @Before
    public void doBefore() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }

    @Test
    public void testAllDupplicatesNoPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ", 0, 1, 100);
            read.setDuplicateReadFlag(true);
            list.add(read);
        }
        Pair<List<SAMRecord>, List<SAMRecord>> myPairing = obj.splitDuplicates(list);
        Assert.assertEquals(0, myPairing.first.size());   // unique
        Assert.assertEquals(10, myPairing.second.size()); // dup's
    }

    @Test
    public void testNoDupplicatesNoPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ", 0, 1, 100);
            read.setDuplicateReadFlag(false);
            list.add(read);
        }
        Pair<List<SAMRecord>, List<SAMRecord>> myPairing = obj.splitDuplicates(list);
        Assert.assertEquals(10, myPairing.first.size()); // unique
        Assert.assertEquals(0, myPairing.second.size()); // dup's
    }

    @Test
    public void testFiftyFiftyNoPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 5; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ", 0, 1, 100);
            read.setDuplicateReadFlag(true);
            list.add(read);
        }
        for (int x = 10; x < 15; x++)
            list.add(ArtificialSAMUtils.createArtificialRead(header, String.valueOf(x), 0, x, 100));

        Pair<List<SAMRecord>, List<SAMRecord>> myPairing = obj.splitDuplicates(list);
        Assert.assertEquals(5, myPairing.first.size());  // unique
        Assert.assertEquals(5, myPairing.second.size()); // dup's
    }

    @Test
    public void testAllDupplicatesAllPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ", 0, 1, 100);
            read.setDuplicateReadFlag(true);
            read.setMateAlignmentStart(100);
            read.setMateReferenceIndex(0);
            read.setReadPairedFlag(true);
            list.add(read);
        }

        Pair<List<SAMRecord>, List<SAMRecord>> myPairing = obj.splitDuplicates(list);
        Assert.assertEquals(0, myPairing.first.size());  // unique
        Assert.assertEquals(10, myPairing.second.size()); // dup's
    }

    @Test
    public void testNoDupplicatesAllPairs() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ", 0, 1, 100);
            if (x == 0) read.setDuplicateReadFlag(true); // one is a dup but (next line)
            read.setMateAlignmentStart(100); // they all have a shared start and mate start so they're dup's
            read.setMateReferenceIndex(0);
            read.setReadPairedFlag(true);
            list.add(read);
        }

        Pair<List<SAMRecord>, List<SAMRecord>> myPairing = obj.splitDuplicates(list);
        Assert.assertEquals(0, myPairing.first.size());  // unique
        Assert.assertEquals(10, myPairing.second.size()); // dup's
    }

    @Test
    public void testAllDupplicatesAllPairsDifferentPairedEnd() {
        List<SAMRecord> list = new ArrayList<SAMRecord>();
        for (int x = 0; x < 10; x++) {
            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "SWEET_READ", 0, 1, 100);
            if (x == 0) read.setDuplicateReadFlag(true); // one is a dup
            read.setMateAlignmentStart(100 + x);
            read.setMateReferenceIndex(0);
            read.setReadPairedFlag(true);
            list.add(read);
        }

        Pair<List<SAMRecord>, List<SAMRecord>> myPairing = obj.splitDuplicates(list);
        Assert.assertEquals(9, myPairing.first.size());  // unique
        Assert.assertEquals(1, myPairing.second.size()); // dup's
    }
}
