package org.broadinstitute.sting.gatk.datasources.providers;

import org.junit.Assert;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import net.sf.samtools.SAMRecord;

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
public class AllLocusViewTest extends LocusViewTemplate {

    @Override
    protected LocusView createView(ShardDataProvider provider) {
        return new AllLocusView(provider);
    }

    /**
     * Test the reads according to an independently derived context.
     * @param view
     * @param range
     * @param reads
     */
    @Override
    protected void testReadsInContext( LocusView view, List<GenomeLoc> range, List<SAMRecord> reads ) {
        AllLocusView allLocusView = (AllLocusView)view;

        // TODO: Should skip over loci not in the given range.
        GenomeLoc firstLoc = range.get(0);
        GenomeLoc lastLoc = range.get(range.size()-1);
        GenomeLoc bounds = GenomeLocParser.createGenomeLoc(firstLoc.getContigIndex(),firstLoc.getStart(),lastLoc.getStop());

        for( long i = bounds.getStart(); i <= bounds.getStop(); i++ ) {
            GenomeLoc site = GenomeLocParser.createGenomeLoc("chr1",i);
            AlignmentContext locusContext = allLocusView.next();
            Assert.assertEquals("Locus context location is incorrect", site, locusContext.getLocation() );
            int expectedReadsAtSite = 0;

            for( SAMRecord read: reads ) {
                if(GenomeLocParser.createGenomeLoc(read).containsP(locusContext.getLocation())) {
                    Assert.assertTrue("Target locus context does not contain reads", locusContext.getReads().contains(read) );
                    expectedReadsAtSite++;
                }
            }

            Assert.assertEquals("Found wrong number of reads at site", expectedReadsAtSite, locusContext.getReads().size());
        }

    }
}
