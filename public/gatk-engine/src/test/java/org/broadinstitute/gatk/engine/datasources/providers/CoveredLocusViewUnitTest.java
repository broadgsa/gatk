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


import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;

import java.util.List;
/**
 * User: hanna
 * Date: May 12, 2009
 * Time: 2:34:46 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Test the CoveredLocusView.
 */
public class CoveredLocusViewUnitTest extends LocusViewTemplate {

    /**
     * Retrieve a covered locus view.
     */
    @Override
    protected LocusView createView(LocusShardDataProvider provider) {
        return new CoveredLocusView(provider);
    }

    /**
     * Test the reads according to an independently derived context.
     * @param view
     * @param range
     * @param reads
     */
    @Override
    protected void testReadsInContext( LocusView view, List<GenomeLoc> range, List<GATKSAMRecord> reads ) {
        CoveredLocusView coveredLocusView = (CoveredLocusView)view;

        // TODO: Should skip over loci not in the given range.
        GenomeLoc firstLoc = range.get(0);
        GenomeLoc lastLoc = range.get(range.size()-1);
        GenomeLoc bounds = genomeLocParser.createGenomeLoc(firstLoc.getContig(),firstLoc.getStart(),lastLoc.getStop());

        for( int i = bounds.getStart(); i <= bounds.getStop(); i++ ) {
            GenomeLoc site = genomeLocParser.createGenomeLoc("chr1",i);

            int expectedReadsAtSite = 0;
            for( GATKSAMRecord read: reads ) {
                if( genomeLocParser.createGenomeLoc(read).containsP(site) )
                    expectedReadsAtSite++;
            }

            if( expectedReadsAtSite < 1 )
                continue;

            Assert.assertTrue(coveredLocusView.hasNext(),"Incorrect number of loci in view");

            AlignmentContext locusContext = coveredLocusView.next();
            Assert.assertEquals(locusContext.getLocation(), site, "Target locus context location is incorrect");
            Assert.assertEquals(locusContext.getReads().size(), expectedReadsAtSite, "Found wrong number of reads at site");

            for( GATKSAMRecord read: reads ) {
                if(genomeLocParser.createGenomeLoc(read).containsP(locusContext.getLocation()))
                    Assert.assertTrue(locusContext.getReads().contains(read),"Target locus context does not contain reads");
            }
        }

        Assert.assertFalse(coveredLocusView.hasNext(),"Iterator is not bounded at boundaries of shard");
    }        
}
