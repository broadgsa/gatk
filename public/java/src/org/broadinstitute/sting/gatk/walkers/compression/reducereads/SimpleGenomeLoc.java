package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.SortedSet;

/**
 * GenomeLocs are very useful objects to keep track of genomic locations and perform set operations
 * with them.
 *
 * However, GenomeLocs are bound to strict validation through the GenomeLocParser and cannot
 * be created easily for small tasks that do not require the rigors of the GenomeLocParser validation
 *
 * SimpleGenomeLoc is a simple utility to create GenomeLocs without going through the parser. Should
 * only be used outside of the engine.
 *
 * User: carneiro
 * Date: 10/16/12
 * Time: 2:07 PM
 */
public class SimpleGenomeLoc extends GenomeLoc {
    private boolean finished;

    public SimpleGenomeLoc(String contigName, int contigIndex, int start, int stop, boolean finished) {
        super(contigName,  contigIndex, start, stop);
        this.finished = finished;
    }

    public boolean isFinished() {
        return finished;
    }

    @Requires("a != null && b != null")
    public static SimpleGenomeLoc merge(SimpleGenomeLoc a, SimpleGenomeLoc b) throws ReviewedStingException {
        if(GenomeLoc.isUnmapped(a) || GenomeLoc.isUnmapped(b)) {
            throw new ReviewedStingException("Tried to merge unmapped genome locs");
        }

        if (!(a.contiguousP(b))) {
            throw new ReviewedStingException("The two genome locs need to be contiguous");
        }


        return new SimpleGenomeLoc(a.getContig(), a.contigIndex,
                Math.min(a.getStart(), b.getStart()),
                Math.max(a.getStop(), b.getStop()),
                a.isFinished());
    }

    /**
     * Merges a list of *sorted* *contiguous* locs into one
     *
     * @param sortedLocs a sorted list of contiguous locs
     * @return one merged loc
     */
    public static SimpleGenomeLoc merge(SortedSet<SimpleGenomeLoc> sortedLocs) {
        SimpleGenomeLoc previousLoc = null;
        for (SimpleGenomeLoc loc : sortedLocs) {
            if (loc.isUnmapped()) {
                throw new ReviewedStingException("Tried to merge unmapped genome locs");
            }
            if (previousLoc != null && !previousLoc.contiguousP(loc)) {
                throw new ReviewedStingException("The genome locs need to be contiguous");
            }
            previousLoc = loc;
        }
        SimpleGenomeLoc firstLoc = sortedLocs.first();
        SimpleGenomeLoc lastLoc = sortedLocs.last();
        return merge(firstLoc, lastLoc);
    }
}
