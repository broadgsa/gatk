package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.LocusWindowWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.*;
import java.io.File;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Apr 23, 2009
 * Time: 10:26:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class TraverseLocusWindows extends TraversalEngine {

    public TraverseLocusWindows(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super(reads, ref, rods);
    }

    public <M,T> T traverse( Walker<M,T> walker,
                             Shard shard,
                             ShardDataProvider dataProvider,
                             T sum ) {

        if ( !(walker instanceof LocusWindowWalker) )
            throw new IllegalArgumentException("Walker isn't a locus window walker!");

        LocusWindowWalker<M, T> locusWindowWalker = (LocusWindowWalker<M, T>)walker;

        GenomeLoc interval = shard.getGenomeLoc();

        ReadView readView = new ReadView( dataProvider );
        LocusReferenceView referenceView = new LocusReferenceView( dataProvider );
        ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );

        LocusContext locus = getLocusContext(readView.iterator(), interval);

        String referenceSubsequence = new String(referenceView.getReferenceBases(interval));

        // Iterate forward to get all reference ordered data covering this interval
        final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation());

        //
        // Execute our contract with the walker.  Call filter, map, and reduce
        //
        final boolean keepMeP = locusWindowWalker.filter(tracker, referenceSubsequence, locus);
        if (keepMeP) {
            M x = locusWindowWalker.map(tracker, referenceSubsequence, locus);
            sum = locusWindowWalker.reduce(x, sum);
        }

        printProgress("intervals", locus.getLocation());

        return sum;
    }

    private LocusContext getLocusContext(Iterator<SAMRecord> readIter, GenomeLoc interval) {
        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        boolean done = false;
        long leftmostIndex = interval.getStart(),
                rightmostIndex = interval.getStop();
        while (readIter.hasNext() && !done) {
            TraversalStatistics.nRecords++;

            SAMRecord read = readIter.next();
            reads.add(read);
            if ( read.getAlignmentStart() < leftmostIndex )
                leftmostIndex = read.getAlignmentStart();
            if ( read.getAlignmentEnd() > rightmostIndex )
                rightmostIndex = read.getAlignmentEnd();
            if ( this.maxReads > 0 && TraversalStatistics.nRecords > this.maxReads ) {
                logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords));
                done = true;
            }
        }

        GenomeLoc window = GenomeLocParser.createGenomeLoc(interval.getContig(), leftmostIndex, rightmostIndex);
        LocusContext locus = new LocusContext(window, reads, null);
        if ( DOWNSAMPLE_BY_COVERAGE )
            locus.downsampleToCoverage(downsamplingCoverage);

        return locus;
    }

}
