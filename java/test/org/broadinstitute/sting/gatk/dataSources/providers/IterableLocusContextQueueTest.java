package org.broadinstitute.sting.gatk.dataSources.providers;

import org.junit.Test;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.LocusShard;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.BaseTest;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileHeader;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Collections;
import java.io.FileNotFoundException;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
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
 * Test the locus context queue.
 */
public class IterableLocusContextQueueTest extends LocusContextQueueTemplate {

    @Override
    protected LocusContextQueue createQueue(ShardDataProvider provider) {
        return new IterableLocusContextQueue(provider);
    }

    /**
     * Test the reads according to an independently derived context.
     * @param queue
     * @param bounds
     * @param reads
     */
    @Override
    protected void testReadsInContext( LocusContextQueue queue, GenomeLoc bounds, List<SAMRecord> reads ) {
        IterableLocusContextQueue iterableQueue = (IterableLocusContextQueue)queue;

        for( long i = bounds.getStart(); i <= bounds.getStop(); i++ ) {
            GenomeLoc site = new GenomeLoc("chr1",i);

            int expectedReadsAtSite = 0;
            for( SAMRecord read: reads ) {
                if( new GenomeLoc(read).containsP(site) )
                    expectedReadsAtSite++;
            }

            if( expectedReadsAtSite < 1 )
                continue;

            Assert.assertTrue("Incorrect number of loci in queue",iterableQueue.hasNext());

            GenomeLoc nextLocus = iterableQueue.next();
            Assert.assertEquals("Next locus context returned is incorrect", site, nextLocus );            

            LocusContext locusContext = iterableQueue.seek(site).peek();
            Assert.assertEquals("Target locus context location is incorrect", site, locusContext.getLocation() );
            Assert.assertEquals("Found wrong number of reads at site", expectedReadsAtSite, locusContext.getReads().size());

            for( SAMRecord read: reads ) {
                if(new GenomeLoc(read).containsP(locusContext.getLocation()))
                    Assert.assertTrue("Target locus context does not contain reads", locusContext.getReads().contains(read) );
            }

        }
    }
}
