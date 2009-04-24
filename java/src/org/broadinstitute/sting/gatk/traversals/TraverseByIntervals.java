package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.IntervalWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;
import java.util.Iterator;
import java.util.ArrayList;
import java.io.File;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import edu.mit.broad.picard.filter.FilteringIterator;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Apr 23, 2009
 * Time: 10:26:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class TraverseByIntervals extends TraversalEngine {

    public TraverseByIntervals(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super(reads, ref, rods);
    }

    public <M,T> T traverse(Walker<M,T> walker, ArrayList<GenomeLoc> locations) {
        if ( walker instanceof IntervalWalker ) {
            IntervalWalker<M, T> intervalWalker = (IntervalWalker<M, T>)walker;
            T sum = traverseByIntervals(intervalWalker, locations);
            return sum;
        } else {
            throw new IllegalArgumentException("Walker isn't an interval walker!");
        }
    }

    /**
     * Traverse by intervals -- the key driver of linearly ordered traversal of intervals.  Provides reads, RODs, and
     * the reference base for each interval in the reference to the intervalWalker walker.  Supports all of the
     * interaction contract implied by the interval walker
     *
     * @param walker An interval walker object
     * @param <M>    MapType -- the result of calling map() on walker
     * @param <T>    ReduceType -- the result of calling reduce() on the walker
     * @return 0 on success
     */
    protected <M, T> T traverseByIntervals(IntervalWalker<M, T> walker, ArrayList<GenomeLoc> locations) {
        logger.debug("Entering traverseByIntervals");

        if(readsFiles.size() > 1)
            throw new UnsupportedOperationException("Cannot do ByInterval traversal on file with multiple inputs.");        

        samReader = initializeSAMFile(readsFiles.get(0));

        verifySortOrder(true);

        walker.initialize();

        T sum = walker.reduceInit();

        if ( locations.isEmpty() ) {
            logger.debug("There are no intervals provided for the traversal");
        } else {
            if ( ! samReader.hasIndex() )
                Utils.scareUser("Processing locations were requested, but no index was found for the input SAM/BAM file. This operation is potentially dangerously slow, aborting.");

            for ( GenomeLoc interval : locations ) {
                logger.debug(String.format("Processing interval %s", interval.toString()));

                CloseableIterator<SAMRecord> readIter = samReader.queryOverlapping( interval.getContig(),
                        (int)interval.getStart(),
                        (int)interval.getStop());

                Iterator<SAMRecord> wrappedIter = WrapReadsIterator( readIter, false );
                sum = carryWalkerOverInterval(walker, wrappedIter, sum, interval);
                readIter.close();
            }
        }

        //printOnTraversalDone("intervals", sum);
        walker.onTraversalDone(sum);
        return sum;
    }

    protected <M, T> T carryWalkerOverInterval(IntervalWalker<M, T> walker, Iterator<SAMRecord> readIter, T sum, GenomeLoc interval ) {
        logger.debug(String.format("TraverseByIntervals.carryWalkerOverInterval Genomic interval is %s", interval));

        // prepare the read filtering read iterator and provide it to a new interval iterator
        FilteringIterator filterIter = new FilteringIterator(readIter, new locusStreamFilterFunc());

        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        ArrayList<Integer> offsets = new ArrayList<Integer>();
        boolean done = false;
        while (filterIter.hasNext() && !done) {
            TraversalStatistics.nRecords++;
            SAMRecord read = filterIter.next();
            reads.add(read);
            offsets.add((int)(read.getAlignmentStart() - interval.getStart()));
            if (this.maxReads > 0 && TraversalStatistics.nRecords > this.maxReads) {
                logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords));
                done = true;
            }
        }

        LocusContext locus = new LocusContext(interval, reads, offsets);
        if ( DOWNSAMPLE_BY_COVERAGE )
            locus.downsampleToCoverage(downsamplingCoverage);

        ReferenceIterator refSite = refIter.seekForward(locus.getLocation());
        locus.setReferenceContig(refSite.getCurrentContig());

        // Iterate forward to get all reference ordered data covering this interval
        final RefMetaDataTracker tracker = getReferenceOrderedDataAtLocus(locus.getLocation());

        sum = walkAtinterval( walker, sum, locus, refSite, tracker );

        //System.out.format("Working at %s\n", locus.getLocation().toString());

        printProgress("intervals", locus.getLocation());

        return sum;
    }

    protected <M, T> T walkAtinterval( final IntervalWalker<M, T> walker,
                                       T sum,
                                       final LocusContext locus,
                                       final ReferenceIterator refSite,
                                       final RefMetaDataTracker tracker ) {
        ReferenceIterator refSiteCopy = refSite;
        StringBuffer refBases = new StringBuffer(refSiteCopy.getBaseAsString());
        int locusLength = (int)(locus.getLocation().getStop() - locus.getLocation().getStart());
        for ( int i = 0; i < locusLength; i++ ) {
            refSiteCopy = refSiteCopy.next();
            refBases.append(refSiteCopy.getBaseAsChar());
        }

        //logger.debug(String.format("  Reference: %s:%d %c", refSite.getCurrentContig().getName(), refSite.getPosition(), refBase));

        //
        // Execute our contract with the walker.  Call filter, map, and reduce
        //
        final boolean keepMeP = walker.filter(tracker, refBases.toString(), locus);
        if (keepMeP) {
            M x = walker.map(tracker, refBases.toString(), locus);
            sum = walker.reduce(x, sum);
        }

        //printProgress("intervals", interval.getLocation());
        return sum;
    }
}