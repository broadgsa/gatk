package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.MockLocusShard;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.features.table.TableCodec;
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Collections;

import net.sf.picard.reference.IndexedFastaSequenceFile;
/**
 * User: hanna
 * Date: May 27, 2009
 * Time: 3:07:23 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Test the transparent view into the reference-ordered data.  At the moment, just do some basic bindings and make
 * sure the data comes through correctly.
 */
public class ReferenceOrderedViewUnitTest extends BaseTest {
    /**
     * Sequence file.
     */
    private static IndexedFastaSequenceFile seq;

    /**
     * our track builder
     */
    RMDTrackBuilder builder = new RMDTrackBuilder();

    @BeforeClass
    public static void init() throws FileNotFoundException {
        // sequence
        seq = new IndexedFastaSequenceFile(new File(hg18Reference));
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    /**
     * Make sure binding to an empty list produces an empty tracker.
     */
    @Test
    public void testNoBindings() {
        Shard shard = new MockLocusShard(Collections.singletonList(GenomeLocParser.createGenomeLoc("chrM",1,30)));
        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, shard.getGenomeLocs().get(0), null, seq, Collections.<ReferenceOrderedDataSource>emptyList());
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(GenomeLocParser.createGenomeLoc("chrM",10));
        Assert.assertEquals("The tracker should not have produced any data", 0, tracker.getAllRods().size());
    }

    /**
     * Test a single ROD binding.
     */
    @Test
    public void testSingleBinding() {
        File file = new File(testDir + "TabularDataTest.dat");
        RMDTrack track = builder.createInstanceOfTrack(TableCodec.class,"tableTest",file);
        ReferenceOrderedDataSource dataSource = new ReferenceOrderedDataSource(null,track);

        Shard shard = new MockLocusShard(Collections.singletonList(GenomeLocParser.createGenomeLoc("chrM",1,30)));

        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, shard.getGenomeLocs().get(0), null, seq, Collections.singletonList(dataSource));
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(GenomeLocParser.createGenomeLoc("chrM",20));
        TableFeature datum = tracker.lookup("tableTest",TableFeature.class);

        Assert.assertEquals("datum parameter for COL1 is incorrect", "C", datum.get("COL1"));
        Assert.assertEquals("datum parameter for COL2 is incorrect", "D", datum.get("COL2"));
        Assert.assertEquals("datum parameter for COL3 is incorrect", "E", datum.get("COL3"));
    }

    /**
     * Make sure multiple bindings are visible from the view.
     */
    @Test
    public void testMultipleBinding() {
        File file = new File(testDir + "TabularDataTest.dat");


        RMDTrack track = builder.createInstanceOfTrack(TableCodec.class,"tableTest1",file);
        ReferenceOrderedDataSource dataSource1 = new ReferenceOrderedDataSource(null,track);
        RMDTrack track2 = builder.createInstanceOfTrack(TableCodec.class,"tableTest2",file);
        ReferenceOrderedDataSource dataSource2 = new ReferenceOrderedDataSource(null,track2);


        Shard shard = new MockLocusShard(Collections.singletonList(GenomeLocParser.createGenomeLoc("chrM",1,30)));

        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, shard.getGenomeLocs().get(0), null, seq, Arrays.asList(dataSource1,dataSource2));
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(GenomeLocParser.createGenomeLoc("chrM",20));
        TableFeature datum1 = tracker.lookup("tableTest1",TableFeature.class);

        Assert.assertEquals("datum1 parameter for COL1 is incorrect", "C", datum1.get("COL1"));
        Assert.assertEquals("datum1 parameter for COL2 is incorrect", "D", datum1.get("COL2"));
        Assert.assertEquals("datum1 parameter for COL3 is incorrect", "E", datum1.get("COL3"));

        TableFeature datum2 = tracker.lookup("tableTest2", TableFeature.class);

        Assert.assertEquals("datum2 parameter for COL1 is incorrect", "C", datum2.get("COL1"));
        Assert.assertEquals("datum2 parameter for COL2 is incorrect", "D", datum2.get("COL2"));
        Assert.assertEquals("datum2 parameter for COL3 is incorrect", "E", datum2.get("COL3"));
    }
}
