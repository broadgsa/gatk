package org.broadinstitute.sting.gatk.dataSources.providers;

import org.junit.Test;
import org.junit.Assert;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.LocusShard;
import org.broadinstitute.sting.gatk.iterators.GenomeLocusIterator;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;
/**
 * User: hanna
 * Date: May 27, 2009
 * Time: 11:10:00 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Tests for viewing the reference from the perspective of a locus.
 */

public class LocusReferenceViewTest extends ReferenceViewTemplate {
    /**
     * Multiple-base pair queries should generate exceptions.
     */
    @Test(expected=InvalidPositionException.class)
    public void testSingleBPFailure() {
        Shard shard = new LocusShard( new GenomeLoc(0,1,50) );

        ShardDataProvider dataProvider = new ShardDataProvider(shard,null,sequenceFile,null);
        LocusReferenceView view = new LocusReferenceView(dataProvider);

        view.getReferenceBase(shard.getGenomeLoc());
    }

    /**
     * Queries outside the bounds of the shard should generate an error.
     */
    @Test(expected=InvalidPositionException.class)
    public void testBoundsFailure() {
        Shard shard = new LocusShard( new GenomeLoc(0,1,50) );

        ShardDataProvider dataProvider = new ShardDataProvider(shard,null,sequenceFile,null);
        LocusReferenceView view = new LocusReferenceView(dataProvider);

        view.getReferenceBase(new GenomeLoc(0,51));
    }


    /**
     * Compares the contents of the fasta and view at a specified location.
     * @param loc
     */
    protected void validateLocation( GenomeLoc loc ) {
        Shard shard = new LocusShard( loc );
        GenomeLocusIterator shardIterator = new GenomeLocusIterator(shard.getGenomeLoc());

        ShardDataProvider dataProvider = new ShardDataProvider(shard,null,sequenceFile,null);
        LocusReferenceView view = new LocusReferenceView(dataProvider);

        while( shardIterator.hasNext() ) {
            GenomeLoc locus = shardIterator.next();

            ReferenceSequence expectedAsSeq = sequenceFile.getSubsequenceAt(locus.getContig(),locus.getStart(),locus.getStop());
            char expected = StringUtil.bytesToString(expectedAsSeq.getBases()).charAt(0);
            char actual = view.getReferenceBase(locus);

            Assert.assertEquals(String.format("Value of base at position %s in shard %s does not match expected",locus.toString(),shard.getGenomeLoc()),
                    expected,
                    actual);
        }
    }

}
