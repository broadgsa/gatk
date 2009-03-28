package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.iterators.LocusIteratorByHanger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;
import java.util.Arrays;
import java.util.Iterator;
import java.util.ArrayList;
import java.io.File;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import edu.mit.broad.picard.filter.FilteringIterator;

/**
 * A simple, short-term solution to iterating over all reference positions over a series of
 * genomic locations. Simply overloads the superclass traverse function to go over the entire
 * interval's reference positions.
 */
public class TraverseByLociByReference extends TraverseByLoci {

    public TraverseByLociByReference(File reads, File ref, List<ReferenceOrderedData> rods) {
        super(reads, ref, rods);
    }

    public <M,T> T traverse(Walker<M,T> walker, ArrayList<GenomeLoc> locations) {
        if ( locations.isEmpty() )
            Utils.scareUser("Requested all locations be processed without providing locations to be processed!");

        return super.traverse(walker, locations);
    }

    protected <M, T> T carryWalkerOverInterval( LocusWalker<M, T> walker,
                                                Iterator<SAMRecord> readIter,
                                                T sum,
                                                GenomeLoc interval ) {
        logger.debug(String.format("TraverseByLociByReference.carryWalkerOverInterval Genomic interval is %s", interval));

        boolean done = false;

        List<SAMRecord> NO_READS = new ArrayList<SAMRecord>();
        List<Integer> NO_OFFSETS = new ArrayList<Integer>();

        FilteringIterator filterIter = new FilteringIterator(readIter, new locusStreamFilterFunc());
        LocusIterator locusIter = new LocusIteratorByHanger(filterIter);        // prepare the iterator by loci from reads
        ReferenceIterator refSite = refIter.seekForward(interval);              // jump to the first reference site
        LocusContext locusFromReads = advanceReadsToLoc(locusIter, interval);    // load up the next locus by reads

        // We keep processing while the next reference location is within the interval
        while ( interval.containsP(refSite.getLocation()) && ! done ) {
            logger.debug(String.format("  LocusFromReads is %s", locusFromReads == null ? null : locusFromReads.getLocation()));

            this.nRecords++;
            GenomeLoc current = refSite.getLocation();
            
            // Iterate forward to get all reference ordered data covering this locus
            final List<ReferenceOrderedDatum> rodData = getReferenceOrderedDataAtLocus(rodIters, current);

            LocusContext locus = null;

            if ( locusFromReads != null && locusFromReads.getLocation().compareTo(current) == 0 ) {   // we are at the same site
                locus = locusFromReads;
                if ( locusIter.hasNext() )
                    locusFromReads = locusIter.next();                          // advance the iterator
            }
            else
                locus = new LocusContext(current, NO_READS, NO_OFFSETS);    // make the empty locus that has no reads

            locus.setReferenceContig(refSite.getCurrentContig());            
            sum = walkAtLocus( walker, sum, locus, refSite, rodData );

            if (this.maxReads > 0 && this.nRecords > this.maxReads) {
                logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + this.nRecords));
                done = true;
            }
            
            refSite = refIter.next();                                       // update our location
        }

        return sum;
    }


    private LocusContext advanceReadsToLoc(LocusIterator locusIter, GenomeLoc target) {
        if ( ! locusIter.hasNext() )
            return null;
        
        LocusContext locus = locusIter.next();

        while ( ! target.containsP(locus.getLocation()) && locusIter.hasNext() ) {
            logger.debug(String.format("  current locus is %s vs %s => %d", locus.getLocation(), target, locus.getLocation().compareTo(target)));
            locus = locusIter.next();
        }

        logger.debug(String.format("  returning %s", locus.getLocation()));
        return locus;
    }
}