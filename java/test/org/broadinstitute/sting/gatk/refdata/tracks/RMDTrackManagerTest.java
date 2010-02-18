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

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.utils.StingException;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

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
            Iterator<GATKFeature> fIter = null;
            fIter = t.getIterator();
            while (fIter.hasNext()) {
                fIter.next();
                count++;
            }
            Assert.assertEquals(100, count);
        }

    }
}

