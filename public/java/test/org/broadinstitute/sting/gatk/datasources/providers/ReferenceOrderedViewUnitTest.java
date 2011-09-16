package org.broadinstitute.sting.gatk.datasources.providers;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.datasources.reads.MockLocusShard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet.RMDStorageType;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

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
    private GenomeLocParser genomeLocParser;

    /**
     * our track builder
     */
    RMDTrackBuilder builder = null;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq);
        builder = new RMDTrackBuilder(seq.getSequenceDictionary(),genomeLocParser,null);
    }

    /**
     * Make sure binding to an empty list produces an empty tracker.
     */
    @Test
    public void testNoBindings() {
        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chrM",1,30)));
        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, genomeLocParser, shard.getGenomeLocs().get(0), null, seq, Collections.<ReferenceOrderedDataSource>emptyList());
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(genomeLocParser.createGenomeLoc("chrM",10), null);
        Assert.assertEquals(tracker.getValues(Feature.class).size(), 0, "The tracker should not have produced any data");
    }

    /**
     * Test a single ROD binding.
     */
    @Test
    public void testSingleBinding() {
        String fileName = testDir + "TabularDataTest.dat";
        RMDTriplet triplet = new RMDTriplet("tableTest","Table",fileName,RMDStorageType.FILE,new Tags());
        ReferenceOrderedDataSource dataSource = new ReferenceOrderedDataSource(triplet,builder,seq.getSequenceDictionary(),genomeLocParser,false);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chrM",1,30)));

        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, genomeLocParser, shard.getGenomeLocs().get(0), null, seq, Collections.singletonList(dataSource));
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(genomeLocParser.createGenomeLoc("chrM",20), null);
        TableFeature datum = tracker.getFirstValue(TableFeature.class, "tableTest");

        Assert.assertEquals(datum.get("COL1"),"C","datum parameter for COL1 is incorrect");
        Assert.assertEquals(datum.get("COL2"),"D","datum parameter for COL2 is incorrect");
        Assert.assertEquals(datum.get("COL3"),"E","datum parameter for COL3 is incorrect");
    }

    /**
     * Make sure multiple bindings are visible from the view.
     */
    @Test
    public void testMultipleBinding() {
        File file = new File(testDir + "TabularDataTest.dat");

        RMDTriplet testTriplet1 = new RMDTriplet("tableTest1","Table",file.getAbsolutePath(),RMDStorageType.FILE,new Tags());
        ReferenceOrderedDataSource dataSource1 = new ReferenceOrderedDataSource(testTriplet1,builder,seq.getSequenceDictionary(),genomeLocParser,false);

        RMDTriplet testTriplet2 = new RMDTriplet("tableTest2","Table",file.getAbsolutePath(),RMDStorageType.FILE,new Tags());
        ReferenceOrderedDataSource dataSource2 = new ReferenceOrderedDataSource(testTriplet2,builder,seq.getSequenceDictionary(),genomeLocParser,false);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chrM",1,30)));

        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, genomeLocParser, shard.getGenomeLocs().get(0), null, seq, Arrays.asList(dataSource1,dataSource2));
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(genomeLocParser.createGenomeLoc("chrM",20), null);
        TableFeature datum1 = tracker.getFirstValue(TableFeature.class, "tableTest1");

        Assert.assertEquals(datum1.get("COL1"),"C","datum1 parameter for COL1 is incorrect");
        Assert.assertEquals(datum1.get("COL2"),"D","datum1 parameter for COL2 is incorrect");
        Assert.assertEquals(datum1.get("COL3"),"E","datum1 parameter for COL3 is incorrect");

        TableFeature datum2 = tracker.getFirstValue(TableFeature.class, "tableTest2");

        Assert.assertEquals(datum2.get("COL1"),"C","datum2 parameter for COL1 is incorrect");
        Assert.assertEquals(datum2.get("COL2"),"D","datum2 parameter for COL2 is incorrect");
        Assert.assertEquals(datum2.get("COL3"),"E","datum2 parameter for COL3 is incorrect");
    }
}
