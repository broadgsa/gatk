package org.broadinstitute.sting.gatk.dataSources.shards;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.dataSources.datum.LocusDatum;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * User: aaron
 * Date: Mar 30, 2009
 * Time: 7:01:56 PM
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
 * @date Mar 30, 2009
 * <p/>
 * Class LociShard
 * <p/>
 * This is the loci shard, which are collectively made when a shatter call is made to
 * a data source.
 */
public class LociShard implements DataShard {

    // our locusIterator
    private final LocusIterator locusIterator;

    // our reference locusIterator
    private final ReferenceIterator refIterator;

    // Iterator over rods
    private final List<ReferenceOrderedData.RODIterator> rodIters;

    // the max number of iterations
    private final int maxCount;

    // how many iterations we've had
    private int iterCount = 0;

    public LociShard(LocusIterator locusIterator, ReferenceIterator refIterator, List<ReferenceOrderedData.RODIterator> rodIters, int maxCount) {
        this.locusIterator = locusIterator;
        this.maxCount = maxCount;
        this.refIterator = refIterator;
        this.rodIters = rodIters;
    }

    public boolean hasNext() {
        return locusIterator.hasNext() && maxCount > iterCount;
    }

    public LocusDatum next() {
        LocusContext locus = locusIterator.next();
        ReferenceIterator refSite = refIterator.seekForward(locus.getLocation());
        locus.setReferenceContig(refSite.getCurrentContig());
        // Iterate forward to get all reference ordered data covering this locus
        final List<ReferenceOrderedDatum> rodData = getReferenceOrderedDataAtLocus(rodIters, locus.getLocation());
        return new LocusDatum(rodData, refSite.getBaseAsChar(), locus);
    }

    public void remove() {
        locusIterator.remove();
    }

    /**
     * Builds a list of the reference ordered datum at loc from each of the iterators.  This function
     * assumes you are accessing the data in order.  You can't use this function for random access.  Each
     * successive call moves you along the file, consuming all data before loc.
     *
     * @param rodIters Iterators to access the RODs
     * @param loc      The location to get the rods at
     * @return A list of ReferenceOrderDatum at loc.  ROD without a datum at loc will be null in the list
     */
    protected List<ReferenceOrderedDatum> getReferenceOrderedDataAtLocus(List<ReferenceOrderedData.RODIterator> rodIters,
                                                                         final GenomeLoc loc) {
        List<ReferenceOrderedDatum> data = new ArrayList<ReferenceOrderedDatum>();
        for (ReferenceOrderedData.RODIterator iter : rodIters) {
            data.add(iter.seekForward(loc));
        }
        return data;
    }
}
