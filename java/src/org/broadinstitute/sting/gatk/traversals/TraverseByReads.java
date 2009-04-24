package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;

import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.LinkedList;
import java.io.File;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 27, 2009
 * Time: 10:26:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class TraverseByReads extends TraversalEngine {

    public TraverseByReads(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super(reads, ref, rods);        
    }

    public <M,T> T traverse(Walker<M,T> walker, ArrayList<GenomeLoc> locations) {
        if ( walker instanceof ReadWalker ) {
            Walker x = walker;
            ReadWalker<?, ?> readWalker = (ReadWalker<?, ?>)x;
            return (T)this.traverseByRead(readWalker, locations);
        } else {
            throw new IllegalArgumentException("Walker isn't a read walker!");
        }
    }

    /**
     * Traverse by read -- the key driver of linearly ordered traversal of reads.  Provides a single read to
     * the walker object, in coordinate order.  Supports all of the
     * interaction contract implied by the read walker
     * sor
     *
     * @param walker A read walker object
     * @param <M>    MapType -- the result of calling map() on walker
     * @param <>    ReduceType -- the result of calling reduce() on the walker
     * @return 0 on success
     */
    public <M, T> Object traverseByRead(ReadWalker<M, T> walker, ArrayList<GenomeLoc> locations) {
        samReadIter = initializeReads();
        if ( refFileName != null && !locations.isEmpty() )
            GenomeLoc.setupRefContigOrdering(new FastaSequenceFile2(refFileName));

        if (refFileName == null && !walker.requiresOrderedReads() && verifyingSamReadIter != null) {
            logger.warn(String.format("STATUS: No reference file provided and unordered reads are tolerated, enabling out of order read processing."));
            if (verifyingSamReadIter != null)
                verifyingSamReadIter.setCheckOrderP(false);
        }

        if ( samReader != null )
            verifySortOrder(refFileName != null || walker.requiresOrderedReads());

        // Initialize the walker
        walker.initialize();

        // Initialize the sum
        T sum = walker.reduceInit();
        List<Integer> offsets = Arrays.asList(0);   // Offset of a single read is always 0

        boolean done = false;
        // copy the locations here in case we ever want to use the full list again later and so that we can remove efficiently
        LinkedList notYetTraversedLocations = new LinkedList(locations);
        while (samReadIter.hasNext() && !done) {
            TraversalStatistics.nRecords++;

            // get the next read
            final SAMRecord read = samReadIter.next();
            final List<SAMRecord> reads = Arrays.asList(read);
            GenomeLoc loc = new GenomeLoc(read);

            // Jump forward in the reference to this locus location
            LocusContext locus = new LocusContext(loc, reads, offsets);
            if (!loc.isUnmapped() && refIter != null) {
                final ReferenceIterator refSite = refIter.seekForward(loc);
                locus.setReferenceContig(refSite.getCurrentContig());
            }

            GenomeLoc.removePastLocs(loc, notYetTraversedLocations);
            if (GenomeLoc.overlapswithSortedLocsP(loc, notYetTraversedLocations, locations.isEmpty())) {

                //
                // execute the walker contact
                //
                final boolean keepMeP = walker.filter(locus, read);
                if (keepMeP) {
                    M x = walker.map(locus, read);
                    sum = walker.reduce(x, sum);
                }

                if (this.maxReads > 0 && TraversalStatistics.nRecords > this.maxReads) {
                    logger.warn(String.format(("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords)));
                    done = true;
                }
            }
            printProgress("reads", loc);

            if (GenomeLoc.pastFinalLocation(loc, locations))
                done = true;
            //System.out.printf("Done? %b%n", done);
        }

        //printOnTraversalDone("reads", sum);
        walker.onTraversalDone(sum);
        return sum;
    }
}
