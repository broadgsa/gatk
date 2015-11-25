/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.engine.datasources.reads.MockLocusShard;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet;
import org.testng.Assert;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.codecs.table.TableFeature;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet.RMDStorageType;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Collections;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
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
        // disable auto-index creation/locking in the RMDTrackBuilder for tests
        builder = new RMDTrackBuilder(seq.getSequenceDictionary(),genomeLocParser,null,true,null);
    }

    /**
     * Make sure binding to an empty list produces an empty tracker.
     */
    @Test
    public void testNoBindings() {
        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chrM",1,30)));
        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, genomeLocParser, shard.getGenomeLocs().get(0), null, seq, Collections.<ReferenceOrderedDataSource>emptyList());
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(genomeLocParser.createGenomeLoc("chrM",10));
        Assert.assertEquals(tracker.getValues(Feature.class).size(), 0, "The tracker should not have produced any data");
    }

    /**
     * Test a single ROD binding.
     */
    @Test
    public void testSingleBinding() {
        String fileName = privateTestDir + "TabularDataTest.dat";
        RMDTriplet triplet = new RMDTriplet("tableTest","Table",fileName,RMDStorageType.FILE,new Tags());
        ReferenceOrderedDataSource dataSource = new ReferenceOrderedDataSource(triplet,builder,seq.getSequenceDictionary(),genomeLocParser,false);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chrM",1,30)));

        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, genomeLocParser, shard.getGenomeLocs().get(0), null, seq, Collections.singletonList(dataSource));
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(genomeLocParser.createGenomeLoc("chrM",20));
        TableFeature datum = tracker.getFirstValue(new RodBinding<TableFeature>(TableFeature.class, "tableTest"));

        Assert.assertEquals(datum.get("COL1"),"C","datum parameter for COL1 is incorrect");
        Assert.assertEquals(datum.get("COL2"),"D","datum parameter for COL2 is incorrect");
        Assert.assertEquals(datum.get("COL3"),"E","datum parameter for COL3 is incorrect");
    }

    /**
     * Make sure multiple bindings are visible from the view.
     */
    @Test
    public void testMultipleBinding() {
        File file = new File(privateTestDir + "TabularDataTest.dat");

        RMDTriplet testTriplet1 = new RMDTriplet("tableTest1","Table",file.getAbsolutePath(),RMDStorageType.FILE,new Tags());
        ReferenceOrderedDataSource dataSource1 = new ReferenceOrderedDataSource(testTriplet1,builder,seq.getSequenceDictionary(),genomeLocParser,false);

        RMDTriplet testTriplet2 = new RMDTriplet("tableTest2","Table",file.getAbsolutePath(),RMDStorageType.FILE,new Tags());
        ReferenceOrderedDataSource dataSource2 = new ReferenceOrderedDataSource(testTriplet2,builder,seq.getSequenceDictionary(),genomeLocParser,false);

        Shard shard = new MockLocusShard(genomeLocParser,Collections.singletonList(genomeLocParser.createGenomeLoc("chrM",1,30)));

        LocusShardDataProvider provider = new LocusShardDataProvider(shard, null, genomeLocParser, shard.getGenomeLocs().get(0), null, seq, Arrays.asList(dataSource1,dataSource2));
        ReferenceOrderedView view = new ManagingReferenceOrderedView( provider );

        RefMetaDataTracker tracker = view.getReferenceOrderedDataAtLocus(genomeLocParser.createGenomeLoc("chrM",20));
        TableFeature datum1 = tracker.getFirstValue(new RodBinding<TableFeature>(TableFeature.class, "tableTest1"));

        Assert.assertEquals(datum1.get("COL1"),"C","datum1 parameter for COL1 is incorrect");
        Assert.assertEquals(datum1.get("COL2"),"D","datum1 parameter for COL2 is incorrect");
        Assert.assertEquals(datum1.get("COL3"),"E","datum1 parameter for COL3 is incorrect");

        TableFeature datum2 = tracker.getFirstValue(new RodBinding<TableFeature>(TableFeature.class, "tableTest2"));

        Assert.assertEquals(datum2.get("COL1"),"C","datum2 parameter for COL1 is incorrect");
        Assert.assertEquals(datum2.get("COL2"),"D","datum2 parameter for COL2 is incorrect");
        Assert.assertEquals(datum2.get("COL3"),"E","datum2 parameter for COL3 is incorrect");
    }
}
