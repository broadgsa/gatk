package org.broadinstitute.sting.gatk.iterators;

import junit.framework.Assert;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * testing of the LocusIteratorByState
 */
public class LocusIteratorByStateTest extends BaseTest {

    private final int MAX_READS = 10;
    private static SAMFileHeader header;
    private LocusIteratorByState li;

    @BeforeClass
    public static void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }


    @Test
    public void testBasicWarnings() {
        // create a bunch of fake records
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        for (int x = 1; x < 50; x++)
            records.add(ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, x, 20));

        // create a test version of the Reads object
        TestReads reads = new TestReads(new ArrayList<File>());
        reads.setMaxPileupSize(MAX_READS);

        // create the iterator by state with the fake reads and fake records
        li = new LocusIteratorByState(records.iterator(), reads);

        // inject the testing version of the locus iterator watcher
        li.setLocusOverflowTracker(new LocusIteratorOverride(MAX_READS));
        ((LocusIteratorOverride)li.getLocusOverflowTracker()).resetWarningCount();
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            //System.err.println(context.getLocation() + " " + context.getPileup().size());
        }
        Assert.assertEquals(1, ((LocusIteratorOverride) li.getLocusOverflowTracker()).getWarningCount());
    }

    @Test
    public void testTwoBigPiles() {
        // create a bunch of fake records
        List<SAMRecord> records = new ArrayList<SAMRecord>();
        for (int x = 1; x < 50; x++)
            records.add(ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, 1, 20));
        for (int x = 1; x < 50; x++)
            records.add(ArtificialSAMUtils.createArtificialRead(header, "readUno", 0, 100, 20));

        // create a test version of the Reads object
        TestReads reads = new TestReads(new ArrayList<File>());
        reads.setMaxPileupSize(MAX_READS);

        // create the iterator by state with the fake reads and fake records
        li = new LocusIteratorByState(records.iterator(), reads);

        // inject the testing version of the locus iterator watcher
        li.setLocusOverflowTracker(new LocusIteratorOverride(MAX_READS));
        ((LocusIteratorOverride)li.getLocusOverflowTracker()).resetWarningCount();
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            //System.err.println(context.getLocation() + " " + context.getPileup().size());
        }
        li.getLocusOverflowTracker().cleanWarningQueue();
        Assert.assertEquals(2, ((LocusIteratorOverride) li.getLocusOverflowTracker()).getWarningCount());
    }
}


class TestReads extends Reads {

    /**
     * Simple constructor for unit testing.
     *
     * @param readsFiles List of reads files to open.
     */
    public TestReads(List<File> readsFiles) {
        super(readsFiles);
    }

    public void setMaxPileupSize(int maxSize) {
        this.maximumReadsAtLocus = maxSize;
    }
}