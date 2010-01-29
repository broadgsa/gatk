package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWindowWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Apr 23, 2009
 * Time: 10:26:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class TraverseLocusWindows extends TraversalEngine {
    /** descriptor of the type */
    private static final String LOCUS_WINDOW_STRING = "intervals";

    public <M,T> T traverse( Walker<M,T> walker,
                             Shard shard,
                             ShardDataProvider dataProvider,
                             T sum ) {

        if ( !(walker instanceof LocusWindowWalker) )
            throw new IllegalArgumentException("Walker isn't a locus window walker!");

        LocusWindowWalker<M, T> locusWindowWalker = (LocusWindowWalker<M, T>)walker;

        if(shard.getGenomeLocs().size() > 1)
            throw new StingException("This traversal does not support multiple intervals within a single shard");
        GenomeLoc interval = shard.getGenomeLocs().get(0);

        ReadView readView = new ReadView( dataProvider );
        LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );
        ReferenceOrderedView referenceOrderedDataView = new ManagingReferenceOrderedView( dataProvider );

         Pair<GenomeLoc, List<SAMRecord>> locus = getLocusContext(readView.iterator(), interval);

        // The TraverseByLocusWindow expands intervals to cover all reads in a non-standard way.
        // TODO: Convert this approach to the standard.
        GenomeLoc expandedInterval = locus.getFirst();

        String referenceSubsequence = new String(referenceView.getReferenceBases(expandedInterval));

        // Iterate forward to get all reference ordered data covering this interval
        final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getFirst());

        //
        // Execute our contract with the walker.  Call filter, map, and reduce
        //
        //final boolean keepMeP = locusWindowWalker.filter(tracker, referenceSubsequence, locus);
        //if (keepMeP) {
        M x = locusWindowWalker.map(tracker, referenceSubsequence, locus.getFirst(), locus.getSecond());
        sum = locusWindowWalker.reduce(x, sum);
        //}

        printProgress(LOCUS_WINDOW_STRING, locus.getFirst());

        return sum;
    }

    private Pair<GenomeLoc, List<SAMRecord>> getLocusContext(StingSAMIterator readIter, GenomeLoc interval) {
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
//        AlignmentContext locus = new AlignmentContext(window, reads, null);
//        if ( readIter.getSourceInfo().getDownsampleToCoverage() != null )
//            locus.downsampleToCoverage(readIter.getSourceInfo().getDownsampleToCoverage());

        return new Pair<GenomeLoc, List<SAMRecord>>(window, reads);
    }

    /**
     * Temporary override of printOnTraversalDone.
     * TODO: Add some sort of TE.getName() function once all TraversalEngines are ported.
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone( T sum ) {
        printOnTraversalDone(LOCUS_WINDOW_STRING, sum );
    }

}
