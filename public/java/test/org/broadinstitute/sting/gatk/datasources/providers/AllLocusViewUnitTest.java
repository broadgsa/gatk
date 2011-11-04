package org.broadinstitute.sting.gatk.datasources.providers;


import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

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
 * Test the view of all loci.
 */
public class AllLocusViewUnitTest extends LocusViewTemplate {

    @Override
    protected LocusView createView(LocusShardDataProvider provider) {
        return new AllLocusView(provider);
    }

    /**
     * Test the reads according to an independently derived context.
     * @param view
     * @param range
     * @param reads
     */
    @Override
    protected void testReadsInContext( LocusView view, List<GenomeLoc> range, List<GATKSAMRecord> reads ) {
        AllLocusView allLocusView = (AllLocusView)view;

        // TODO: Should skip over loci not in the given range.
        GenomeLoc firstLoc = range.get(0);
        GenomeLoc lastLoc = range.get(range.size()-1);
        GenomeLoc bounds = genomeLocParser.createGenomeLoc(firstLoc.getContig(),firstLoc.getStart(),lastLoc.getStop());

        for( int i = bounds.getStart(); i <= bounds.getStop(); i++ ) {
            GenomeLoc site = genomeLocParser.createGenomeLoc("chr1",i);
            AlignmentContext locusContext = allLocusView.next();
            Assert.assertEquals(locusContext.getLocation(), site, "Locus context location is incorrect");
            int expectedReadsAtSite = 0;

            for( GATKSAMRecord read: reads ) {
                if(genomeLocParser.createGenomeLoc(read).containsP(locusContext.getLocation())) {
                    Assert.assertTrue(locusContext.getReads().contains(read),"Target locus context does not contain reads");
                    expectedReadsAtSite++;
                }
            }

            Assert.assertEquals(locusContext.getReads().size(), expectedReadsAtSite, "Found wrong number of reads at site");
        }

    }
}
