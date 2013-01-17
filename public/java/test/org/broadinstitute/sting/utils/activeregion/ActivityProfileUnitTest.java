/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.activeregion;


// the imports for unit testing.


import net.sf.picard.reference.ReferenceSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


public class ActivityProfileUnitTest extends BaseTest {
    private GenomeLocParser genomeLocParser;
    private GenomeLoc startLoc;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq);
        startLoc = genomeLocParser.createGenomeLoc("chr1", 1, 1, 100);
    }

    // --------------------------------------------------------------------------------
    //
    // Basic tests Provider
    //
    // --------------------------------------------------------------------------------

    private class BasicActivityProfileTestProvider extends TestDataProvider {
        List<Double> probs;
        List<ActiveRegion> expectedRegions;
        int extension = 0;
        GenomeLoc regionStart = startLoc;
        final ProfileType type;

        public BasicActivityProfileTestProvider(final ProfileType type, final List<Double> probs, boolean startActive, int ... startsAndStops) {
            super(BasicActivityProfileTestProvider.class);
            this.type = type;
            this.probs = probs;
            this.expectedRegions = toRegions(startActive, startsAndStops);
            setName(getName());
        }

        private String getName() {
            return String.format("type=%s probs=%s expectedRegions=%s", type, Utils.join(",", probs), Utils.join(",", expectedRegions));
        }

        public ActivityProfile makeProfile() {
            switch ( type ) {
                case Base: return new ActivityProfile(genomeLocParser);
                case BandPass: return new BandPassActivityProfile(genomeLocParser);
                default: throw new IllegalStateException(type.toString());
            }
        }

        private List<ActiveRegion> toRegions(boolean isActive, int[] startsAndStops) {
            List<ActiveRegion> l = new ArrayList<ActiveRegion>();
            for ( int i = 0; i < startsAndStops.length - 1; i++) {
                int start = regionStart.getStart() + startsAndStops[i];
                int end = regionStart.getStart() + startsAndStops[i+1] - 1;
                GenomeLoc activeLoc = genomeLocParser.createGenomeLoc(regionStart.getContig(), start, end);
                ActiveRegion r = new ActiveRegion(activeLoc, isActive, genomeLocParser, extension);
                l.add(r);
                isActive = ! isActive;
            }
            return l;
        }
    }

    private enum ProfileType {
        Base, BandPass
    }

    @DataProvider(name = "BasicActivityProfileTestProvider")
    public Object[][] makeQualIntervalTestProvider() {
        for ( final ProfileType type : ProfileType.values() ) {
            new BasicActivityProfileTestProvider(type, Arrays.asList(1.0), true, 0, 1);
            new BasicActivityProfileTestProvider(type, Arrays.asList(1.0, 0.0), true, 0, 1, 2);
            new BasicActivityProfileTestProvider(type, Arrays.asList(0.0, 1.0), false, 0, 1, 2);
            new BasicActivityProfileTestProvider(type, Arrays.asList(1.0, 0.0, 1.0), true, 0, 1, 2, 3);
            new BasicActivityProfileTestProvider(type, Arrays.asList(1.0, 1.0, 1.0), true, 0, 3);
        }

        return BasicActivityProfileTestProvider.getTests(BasicActivityProfileTestProvider.class);
    }

    @Test(dataProvider = "BasicActivityProfileTestProvider")
    public void testBasicActivityProfile(BasicActivityProfileTestProvider cfg) {
        ActivityProfile profile = cfg.makeProfile();

        Assert.assertTrue(profile.isEmpty());

        Assert.assertEquals(profile.parser, genomeLocParser);

        for ( int i = 0; i < cfg.probs.size(); i++ ) {
            double p = cfg.probs.get(i);
            GenomeLoc loc = genomeLocParser.createGenomeLoc(cfg.regionStart.getContig(), cfg.regionStart.getStart() + i, cfg.regionStart.getStart() + i);
            profile.add(new ActivityProfileState(loc, p));
            Assert.assertFalse(profile.isEmpty());
        }
        Assert.assertEquals(profile.regionStartLoc, genomeLocParser.createGenomeLoc(cfg.regionStart.getContig(), cfg.regionStart.getStart(), cfg.regionStart.getStart() ));

        Assert.assertEquals(profile.size(), cfg.probs.size());
        assertProbsAreEqual(profile.isActiveList, cfg.probs);

        assertRegionsAreEqual(profile.createActiveRegions(0, 100), cfg.expectedRegions);

        Assert.assertEquals(profile.createDerivedProfile(profile.isActiveList).getClass(), profile.getClass());

        final List<ActivityProfileState> empty = new LinkedList<ActivityProfileState>();
        Assert.assertEquals(profile.createDerivedProfile(empty).size(), 0);
    }

    private void assertRegionsAreEqual(List<ActiveRegion> actual, List<ActiveRegion> expected) {
        Assert.assertEquals(actual.size(), expected.size());
        for ( int i = 0; i < actual.size(); i++ ) {
            Assert.assertTrue(actual.get(i).equalExceptReads(expected.get(i)));
        }
    }

    private void assertProbsAreEqual(List<ActivityProfileState> actual, List<Double> expected) {
        Assert.assertEquals(actual.size(), expected.size());
        for ( int i = 0; i < actual.size(); i++ ) {
            Assert.assertEquals(actual.get(i).isActiveProb, expected.get(i));
        }
    }

    // todo -- test extensions
}