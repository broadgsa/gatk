package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.providers.LocusReferenceView;
import org.broadinstitute.sting.gatk.datasources.providers.ReadView;
import org.broadinstitute.sting.gatk.datasources.providers.ReferenceOrderedView;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWindowWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Apr 23, 2009
 * Time: 10:26:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class TraverseLocusWindows extends TraversalEngine {

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

        AlignmentContext locus = getLocusContext(readView.iterator(), interval);

        // The TraverseByLocusWindow expands intervals to cover all reads in a non-standard way.
        // TODO: Convert this approach to the standard.
        GenomeLoc expandedInterval = locus.getLocation();

        String referenceSubsequence = new String(referenceView.getReferenceBases(expandedInterval));

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

    private AlignmentContext getLocusContext(StingSAMIterator readIter, GenomeLoc interval) {
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
            if ( this.maximumIterations > 0 && TraversalStatistics.nRecords > this.maximumIterations) {
                logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords));
                done = true;
            }
        }

        GenomeLoc window = GenomeLocParser.createGenomeLoc(interval.getContig(), leftmostIndex, rightmostIndex);
        AlignmentContext locus = new AlignmentContext(window, reads, null);
        if ( readIter.getSourceInfo().getDownsampleToCoverage() != null )
            locus.downsampleToCoverage(readIter.getSourceInfo().getDownsampleToCoverage());

        return locus;
    }

    /**
     * Temporary override of printOnTraversalDone.
     * TODO: Add some sort of TE.getName() function once all TraversalEngines are ported.
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone( T sum ) {
        printOnTraversalDone( "intervals", sum );
    }

}
