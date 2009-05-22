package org.broadinstitute.sting.gatk.dataSources.providers;

import org.junit.Assert;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.LocusContext;
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
 * Test the CoveredLocusView.
 */
public class CoveredLocusViewTest extends LocusViewTemplate {

    /**
     * Retrieve a covered locus view.
     */
    @Override
    protected LocusView createView(ShardDataProvider provider) {
        return new CoveredLocusView(provider);
    }

    /**
     * Test the reads according to an independently derived context.
     * @param view
     * @param bounds
     * @param reads
     */
    @Override
    protected void testReadsInContext( LocusView view, GenomeLoc bounds, List<SAMRecord> reads ) {
        CoveredLocusView coveredLocusView = (CoveredLocusView)view;

        for( long i = bounds.getStart(); i <= bounds.getStop(); i++ ) {
            GenomeLoc site = new GenomeLoc("chr1",i);

            int expectedReadsAtSite = 0;
            for( SAMRecord read: reads ) {
                if( new GenomeLoc(read).containsP(site) )
                    expectedReadsAtSite++;
            }

            if( expectedReadsAtSite < 1 )
                continue;

            Assert.assertTrue("Incorrect number of loci in view",coveredLocusView.hasNext());

            LocusContext locusContext = coveredLocusView.next();
            Assert.assertEquals("Target locus context location is incorrect", site, locusContext.getLocation() );
            Assert.assertEquals("Found wrong number of reads at site", expectedReadsAtSite, locusContext.getReads().size());

            for( SAMRecord read: reads ) {
                if(new GenomeLoc(read).containsP(locusContext.getLocation()))
                    Assert.assertTrue("Target locus context does not contain reads", locusContext.getReads().contains(read) );
            }
        }

        Assert.assertFalse("Iterator is not bounded at boundaries of shard", coveredLocusView.hasNext());
    }        
}
