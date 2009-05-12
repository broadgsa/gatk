package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;
import java.util.ArrayList;
import java.io.File;

import net.sf.samtools.SAMRecord;

/**
 * A simple, short-term solution to iterating over all reference positions over a series of
 * genomic locations. Simply overloads the superclass traverse function to go over the entire
 * interval's reference positions.
 */
public class TraverseByReference extends TraverseByLoci {
    
    public TraverseByReference(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super(reads, ref, rods);

        logger.debug("Creating TraverseByReference");

        if ( reads != null )
            throw new IllegalArgumentException("By reference traversal doesn't support reads, but file was given " + reads);
        if ( ref == null )
            throw new IllegalArgumentException("By reference traversal requires reference file but none was given");
    }

    /**
     * Traverse by loci -- the key driver of linearly ordered traversal of loci.  Provides reads, RODs, and
     * the reference base for each locus in the reference to the LocusWalker walker.  Supports all of the
     * interaction contract implied by the locus walker
     *
     * @param walker A locus walker object
     * @param <M>    MapType -- the result of calling map() on walker
     * @param <T>    ReduceType -- the result of calling reduce() on the walker
     * @return 0 on success
     */
    protected <M, T> T traverseByLoci(LocusWalker<M, T> walker, List<GenomeLoc> locations) {
        logger.debug("Entering traverseByReference");

        // initialize the walker object
        walker.initialize();

        T sum = walker.reduceInit();
        if ( ! locations.isEmpty() ) {
            logger.debug("Doing interval-based traversal");

            // we are doing interval-based traversals
            for ( GenomeLoc interval : locations ) {
                logger.debug(String.format("Processing locus %s", interval.toString()));
                sum = carryWalkerOverReference(walker, sum, interval);
            }
        }
        else {
            // We aren't locus oriented
            logger.debug("Doing non-interval-based traversal");
            sum = carryWalkerOverReference(walker, sum, null);
        }

        //printOnTraversalDone("reference", sum);
        walker.onTraversalDone(sum);
        return sum;
    }

    protected <M, T> T carryWalkerOverReference( LocusWalker<M, T> walker,
                                                 T sum,
                                                 GenomeLoc interval ) {
        logger.debug(String.format("traverseByReference.carryWalkerOverReference Genomic interval is %s", interval));

        boolean done = false;

        List<SAMRecord> NO_READS = new ArrayList<SAMRecord>();
        List<Integer> NO_OFFSETS = new ArrayList<Integer>();

        ReferenceIterator refSite = null;
        if ( interval != null )
            refSite = refIter.seekForward(interval);              // jump to the first reference site
        else
            refSite = refIter.next();
        
        // We keep processing while the next reference location is within the interval
        while ( (interval == null || interval.containsP(refSite.getLocation())) && ! done ) {
            TraversalStatistics.nRecords++;
            GenomeLoc current = refSite.getLocation();

            // Iterate forward to get all reference ordered data covering this locus
            final RefMetaDataTracker tracker = getReferenceOrderedDataAtLocus(current);

            LocusContext locus = new LocusContext(current, NO_READS, NO_OFFSETS);    // make the empty locus that has no reads
            locus.setReferenceContig(refSite.getCurrentContig());
            sum = walkAtLocus( walker, sum, locus, refSite, tracker );

            if (this.maxReads > 0 && TraversalStatistics.nRecords > this.maxReads) {
                logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords));
                done = true;
            }

            //printProgress("ref", locus.getLocation());
            refSite = refIter.next();                                       // update our location
        }

        return sum;
    }
}