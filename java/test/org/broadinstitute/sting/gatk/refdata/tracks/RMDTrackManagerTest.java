/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.tracks;

import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * class RMDTrackManagerTest
 * tests out the ability of the RMDTrackManager to correctly create RMDtracks based on the requested types.
 */
public class RMDTrackManagerTest extends BaseTest {
    List<String> triplets;
    List<RMDTrack> tracks;

    @Before
    public void setup() {
        RMDTrackManager manager = new RMDTrackManager();
        triplets = new ArrayList<String>();

        // add our db snp data
        triplets.add("MyDbSNP");
        triplets.add("DBSNP");
        triplets.add("testdata/small.dbsnp.rod");
        tracks = manager.getReferenceMetaDataSources(triplets);
    }

    @Test
    public void testBuilderQuery() {
        for (RMDTrack t : tracks) {
            System.err.println("name = " + t.getName() + " type = " + t.getType().getSimpleName() + " file = " + t.getFile());
            int count = 0;
            Iterator<GATKFeature> fIter;
            try {
                fIter = ((FeatureReaderTrack) t).query("1", 1, 5000);
            } catch (IOException e) {
                throw new StingException("blah I/O exception");
            }
            while (fIter.hasNext()) {
                fIter.next();
                count++;
            }
            Assert.assertEquals(100, count);
        }
    }

    @Test
    public void testBuilderIterator() {
        for (RMDTrack t : tracks) {
            System.err.println("name = " + t.getName() + " type = " + t.getType().getSimpleName() + " file = " + t.getFile());
            int count = 0;
            Iterator<GATKFeature> fIter = t.getIterator();
            while (fIter.hasNext()) {
                fIter.next();
                count++;
            }
            Assert.assertEquals(100, count);
        }
    }

    // @Test used only to determine how fast queries are, don't uncomment! (unless you know what you're doing).
    public void testSpeedOfRealQuery() {
        IndexedFastaSequenceFile file = null;
        try {
            file = new IndexedFastaSequenceFile(new File("/broad/1KG/reference/human_b36_both.fasta"));
        } catch (FileNotFoundException e) {
            Assert.assertTrue(false);
        }
        final int intervalSize = 10000000;
        GenomeLocParser.setupRefContigOrdering(file.getSequenceDictionary());
        RMDTrackManager manager = new RMDTrackManager();
        // add our db snp data
        triplets.clear();
        triplets.add("db");
        triplets.add("DBSNP");
        triplets.add("../../GATK_Data/dbsnp_130_b36.rod");
        Assert.assertEquals(1, manager.getReferenceMetaDataSources(triplets).size());
        RMDTrack t = manager.getReferenceMetaDataSources(triplets).get(0);
        // make sure we have a single track
        // lets test the first and 20th contigs of the human reference

        for (int loop = 1; loop <= 22; loop++) {
            SAMSequenceRecord seqRec = GenomeLocParser.getContigInfo(String.valueOf(loop));
            String name = seqRec.getSequenceName();
            Iterator<GATKFeature> fIter;
            for (int x = 1; x < seqRec.getSequenceLength() - intervalSize; x += intervalSize) {
                long firstTime = System.currentTimeMillis();
                long count = 0;
                try {
                    fIter = ((FeatureReaderTrack) t).query("1", x, x + intervalSize);
                } catch (IOException e) {
                    throw new StingException("blah I/O exception");
                }
                while (fIter.hasNext()) {
                    fIter.next();
                    count++;
                }
                System.err.println(name + "," + count + "," + (System.currentTimeMillis() - firstTime));
            }
        }
    }
}

