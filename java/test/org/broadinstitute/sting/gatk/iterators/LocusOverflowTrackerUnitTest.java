package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;


/**
 * @author aaron
 *         <p/>
 *         Class LocusOverflowTrackerUnitTest
 *         <p/>
 *         test out the locus overflow tracker
 */
public class LocusOverflowTrackerUnitTest extends BaseTest {

    private LocusOverflowTracker tracker;
    private final int MAX_READS = 10;
    private static SAMFileHeader header;

    @BeforeClass
    public static void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }

    @Before
    public void beforeTest() {
        tracker = new LocusIteratorOverride(MAX_READS);
        ((LocusIteratorOverride) tracker).resetWarningCount();
    }

    @Test
    public void testLocusOverflow() {
        SAMRecord rec = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, 1, 100);
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(rec);
        if (tracker.exceeded(loc, MAX_READS - 1))
            Assert.fail("We shouldn't be exceeded when MAX_READS -1 is the input");
        if (!tracker.exceeded(loc, MAX_READS)) Assert.fail("We should be exceeded when MAX_READS is the input");
        if (!tracker.exceeded(loc, MAX_READS + 1))
            Assert.fail("We shouldn't be exceeded when MAX_READS +1 is the input");
    }

    @Test
    public void testContinuousLocus() {
        for (int x = 1; x < 5; x++) {
            SAMRecord rec = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, x, 100);
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(rec);
            tracker.exceeded(loc, MAX_READS + 1);
        }
        tracker.cleanWarningQueue();
        Assert.assertEquals(1, ((LocusIteratorOverride) tracker).getWarningCount());
    }

    @Test
    public void testTwoSeperateContinuousLoci() {
        for (int x = 1; x < 5; x++) {
            SAMRecord rec = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, x, 2);
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(rec);
            tracker.exceeded(loc, MAX_READS + 1);
        }
        for (int x = 10; x < 15; x++) {
            SAMRecord rec = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, x, 2);
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(rec);
            tracker.exceeded(loc, MAX_READS + 1);
        }
        tracker.cleanWarningQueue();
        Assert.assertEquals(2, ((LocusIteratorOverride) tracker).getWarningCount());
    }

    @Test
    // make sure we get only the specified number of warnings
    public void testOverflow() {
        for (int x = 1; x < (LocusOverflowTracker.warningsEmitted * 3); x += 2) {
            SAMRecord rec = ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, x, 100);
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(rec);
            tracker.exceeded(loc, MAX_READS + 1);
        }
        tracker.cleanWarningQueue();
        Assert.assertEquals(LocusOverflowTracker.warningsEmitted, ((LocusIteratorOverride) tracker).getWarningCount());
    }

}


class LocusIteratorOverride extends LocusOverflowTracker {

    public int warningsDropped = 0;

    public LocusIteratorOverride(int maxPileup) {
        super(maxPileup);
    }

    /** override to count warnings instead of generating output. */
    protected void warnUser() {
        warningInQueue = false;
        if (warningsEmitted <= MAX_WARNINGS)
            warningsEmitted++;
        else
            warningsDropped++;
    }

    public int getWarningCount() {
        return warningsEmitted;
    }

    public void resetWarningCount() {
        LocusOverflowTracker.warningsEmitted = 0;
    }
}