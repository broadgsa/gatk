package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;

/**
 *
 * User: aaron
 * Date: Apr 2, 2009
 * Time: 2:12:12 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 2, 2009
 * <p/>
 * Class BoundedReferenceIterator
 * <p/>
 * This class is a decorator class from Reference Iterator (though it is constrained
 * by the fact that referenceIterator.seekForwardOffset explicitly returns a referenceIterator
 * for now
 * <p/>
 * TODO: Fix the underlying iterator and this class to model a real decorator pattern
 */
public class BoundedReferenceIterator implements Iterator<ReferenceIterator> {
    // the location to screen over
    private final GenomeLoc mLoc;
    private final ReferenceIterator referenceIterator;

    /**
     * Default constructor
     *
     * @param referenceIterator
     * @param loc
     */
    public BoundedReferenceIterator(ReferenceIterator referenceIterator, GenomeLoc loc) {
        this.referenceIterator = referenceIterator;
        this.mLoc = loc;
    }


    /**
     * isSubRegion
     * <p/>
     * returns true if we include the whole passed in region
     *
     * @param loc the genome region to check
     * @return true if we include THE WHOLE specified region
     */
    protected boolean isSubRegion(GenomeLoc loc) {
        // if the location is null, we assume we're all inclusive (we represent the whole genome).
        if (mLoc == null || loc.isBetween(mLoc, mLoc)) {
            return true;
        }
        return false;
    }

    /**
     * returns true if we include the whole passed in region
     *
     * @param contig
     * @param start
     * @param stop
     * @return true if we enclose the passed region, false otherwise
     */
    protected boolean isSubRegion(final String contig, final int start, final int stop) {
        final GenomeLoc lc = new GenomeLoc(contig, start, stop);
        return isSubRegion(lc);
    }

    /**
     * If we're less then the limiting genomeLoc
     *
     * @param loc
     * @return
     */
    protected boolean isLessThan(GenomeLoc loc) {
        return loc.isPast(mLoc);
    }


    // our adapted next function
    public boolean hasNext() {
        // first check that we are within the search place
        GenomeLoc loc = referenceIterator.getLocation();
        if (!isSubRegion(loc)) {
            return false;
        }

        return referenceIterator.hasNext();
    }

    public ReferenceIterator next() {
        return referenceIterator.next();
    }

    public void remove() {
        referenceIterator.remove();
    }


}
