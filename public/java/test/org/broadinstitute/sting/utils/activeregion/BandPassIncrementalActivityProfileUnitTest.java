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
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


public class BandPassIncrementalActivityProfileUnitTest extends BaseTest {
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    @DataProvider(name = "BandPassBasicTest")
    public Object[][] makeBandPassTest() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        for ( int start : Arrays.asList(1, 10, 100, 1000) ) {
            for ( boolean precedingIsActive : Arrays.asList(true, false) ) {
                for ( int precedingSites: Arrays.asList(0, 1, 10, 100) ) {
                    for ( int bandPassSize : Arrays.asList(0, 1, 10, 100) ) {
//        for ( int start : Arrays.asList(10) ) {
//            for ( boolean precedingIsActive : Arrays.asList(false) ) {
//                for ( int precedingSites: Arrays.asList(0) ) {
//                    for ( int bandPassSize : Arrays.asList(1) ) {
                        tests.add(new Object[]{ start, precedingIsActive, precedingSites, bandPassSize });
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BandPassBasicTest")
    public void testBandPass(final int start, final boolean precedingIsActive, final int nPrecedingSites, final int bandPassSize) {
        final BandPassIncrementalActivityProfile profile = new BandPassIncrementalActivityProfile(genomeLocParser, bandPassSize);

        final int expectedBandSize = bandPassSize * 2 + 1;
        Assert.assertEquals(profile.getBandSize(), expectedBandSize, "Wrong expected band size");

        final String contig = genomeLocParser.getContigs().getSequences().get(0).getSequenceName();
        final double precedingProb = precedingIsActive ? 1.0 : 0.0;
        for ( int i = 0; i < nPrecedingSites; i++ ) {
            final GenomeLoc loc = genomeLocParser.createGenomeLoc(contig, i + start);
            final ActivityProfileState state = new ActivityProfileState(loc, precedingProb);
            profile.add(state);
        }

        final GenomeLoc nextLoc = genomeLocParser.createGenomeLoc(contig, nPrecedingSites + start);
        profile.add(new ActivityProfileState(nextLoc, 1.0));

        if ( precedingIsActive == false && nPrecedingSites >= bandPassSize && bandPassSize < start ) {
            // we have enough space that all probs fall on the genome
            final double[] probs = profile.getProbabilitiesAsArray();
            Assert.assertEquals(MathUtils.sum(probs), 1.0 * (nPrecedingSites * precedingProb + 1), 1e-3, "Activity profile doesn't sum to number of non-zero prob states");
        }
    }

    private double[] bandPassInOnePass(final BandPassIncrementalActivityProfile profile, final double[] activeProbArray) {
        final double[] bandPassProbArray = new double[activeProbArray.length];

        // apply the band pass filter for activeProbArray into filteredProbArray
        final double[] GaussianKernel = profile.getKernel();
        for( int iii = 0; iii < activeProbArray.length; iii++ ) {
            final double[] kernel = ArrayUtils.subarray(GaussianKernel, Math.max(profile.getFilteredSize() - iii, 0), Math.min(GaussianKernel.length, profile.getFilteredSize() + activeProbArray.length - iii));
            final double[] activeProbSubArray = ArrayUtils.subarray(activeProbArray, Math.max(0,iii - profile.getFilteredSize()), Math.min(activeProbArray.length,iii + profile.getFilteredSize() + 1));
            bandPassProbArray[iii] = MathUtils.dotProduct(activeProbSubArray, kernel);
        }

        return bandPassProbArray;
    }

    @DataProvider(name = "BandPassComposition")
    public Object[][] makeBandPassComposition() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        for ( int bandPassSize : Arrays.asList(0, 1, 10, 100, BandPassIncrementalActivityProfile.DEFAULT_FILTER_SIZE) ) {
            for ( int integrationLength : Arrays.asList(1, 10, 100, 1000) ) {
                tests.add(new Object[]{ bandPassSize, integrationLength });
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test( dataProvider = "BandPassComposition")
    public void testBandPassComposition(final int bandPassSize, final int integrationLength) {
        final int start = 1;
        final BandPassIncrementalActivityProfile profile = new BandPassIncrementalActivityProfile(genomeLocParser, bandPassSize);
        final double[] rawActiveProbs = new double[integrationLength + bandPassSize * 2];

        // add a buffer so that we can get all of the band pass values
        final String contig = genomeLocParser.getContigs().getSequences().get(0).getSequenceName();
        int pos = start;
        int rawProbsOffset = 0;
        for ( int i = 0; i < bandPassSize; i++ ) {
            final GenomeLoc loc = genomeLocParser.createGenomeLoc(contig, pos++);
            final ActivityProfileState state = new ActivityProfileState(loc, 0.0);
            profile.add(state);
            rawActiveProbs[rawProbsOffset++] = 0.0;
            rawActiveProbs[rawActiveProbs.length - rawProbsOffset] = 0.0;
        }

        for ( int i = 0; i < integrationLength; i++ ) {
            final GenomeLoc nextLoc = genomeLocParser.createGenomeLoc(contig, pos++);
            profile.add(new ActivityProfileState(nextLoc, 1.0));
            rawActiveProbs[rawProbsOffset++] = 1.0;

            for ( int j = 0; j < profile.size(); j++ ) {
                Assert.assertTrue(profile.getStateList().get(j).isActiveProb >= 0.0, "State probability < 0 at " + j);
                Assert.assertTrue(profile.getStateList().get(j).isActiveProb <= 1.0 + 1e-3, "State probability > 1 at " + j);
            }
        }

        final double[] expectedProbs = bandPassInOnePass(profile, rawActiveProbs);
        for ( int j = 0; j < profile.size(); j++ ) {
            Assert.assertEquals(profile.getStateList().get(j).isActiveProb, expectedProbs[j], "State probability not expected at " + j);
        }
    }
}
