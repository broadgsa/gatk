package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.LocusReferenceView;
import org.broadinstitute.sting.gatk.datasources.providers.LocusView;
import org.broadinstitute.sting.gatk.datasources.providers.ReferenceOrderedView;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * A simple solution to iterating over all reference positions over a series of genomic locations.
 */
public class TraverseLociLinear<M,T> extends TraverseLociBase<M,T> {

    @Override
    protected TraverseResults<T> traverse(LocusWalker<M, T> walker, LocusView locusView, LocusReferenceView referenceView, ReferenceOrderedView referenceOrderedDataView, T sum) {
        // We keep processing while the next reference location is within the interval
        boolean done = false;
        int numIterations = 0;

        while( locusView.hasNext() && ! done ) {
            numIterations++;
            final AlignmentContext locus = locusView.next();
            final GenomeLoc location = locus.getLocation();

            // create reference context. Note that if we have a pileup of "extended events", the context will
            // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
            final ReferenceContext refContext = referenceView.getReferenceContext(location);

            // Iterate forward to get all reference ordered data covering this location
            final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation(), refContext);

            final boolean keepMeP = walker.filter(tracker, refContext, locus);
            if (keepMeP) {
                final M x = walker.map(tracker, refContext, locus);
                sum = walker.reduce(x, sum);
                done = walker.isDone();
            }

            printProgress(locus.getLocation());
        }

        return new TraverseResults<T>(numIterations, sum);
    }
}
