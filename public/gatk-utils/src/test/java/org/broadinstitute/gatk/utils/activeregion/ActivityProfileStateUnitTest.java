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

package org.broadinstitute.gatk.utils.activeregion;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 1/17/13
 * Time: 2:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class ActivityProfileStateUnitTest {
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 100);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @DataProvider(name = "ActiveProfileResultProvider")
    public Object[][] makeActiveProfileResultProvider() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        final String chr = genomeLocParser.getContigs().getSequence(0).getSequenceName();
        for ( final GenomeLoc loc : Arrays.asList(
                genomeLocParser.createGenomeLoc(chr, 10, 10),
                genomeLocParser.createGenomeLoc(chr, 100, 100) )) {
            for ( final double prob : Arrays.asList(0.0, 0.5, 1.0) ) {
                for ( final ActivityProfileState.Type state : ActivityProfileState.Type.values() ) {
                    for ( final Number value : Arrays.asList(1, 2, 4) ) {
                        tests.add(new Object[]{ loc, prob, state, value});
                    }
                }
                tests.add(new Object[]{ loc, prob, null, null});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ActiveProfileResultProvider")
    public void testActiveProfileResultProvider(GenomeLoc loc, final double prob, ActivityProfileState.Type maybeState, final Number maybeNumber) {
        final ActivityProfileState apr = maybeState == null
                ? new ActivityProfileState(loc, prob)
                : new ActivityProfileState(loc, prob, maybeState, maybeNumber);

        Assert.assertEquals(apr.getLoc(), loc);
        Assert.assertNotNull(apr.toString());
        Assert.assertEquals(apr.isActiveProb, prob);
        Assert.assertEquals(apr.resultState, maybeState == null ? ActivityProfileState.Type.NONE : maybeState);
        Assert.assertEquals(apr.resultValue, maybeState == null ? null : maybeNumber);
    }
}
