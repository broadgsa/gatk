package org.broadinstitute.sting.gatk.dataSources.datum;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 *
 * User: aaron
 * Date: Mar 30, 2009
 * Time: 3:08:28 PM
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
 * Class LocusDatum
 * <p/>
 * The datum for loci.  It contains the reference base, locusContext,
 * and the reference order data.
 */
public class LocusDatum implements Datum {

    // our reference order data    
    private final List<ReferenceOrderedDatum> rodData;
    // our seq base
    private final char ref;
    // our locus context  
    private final LocusContext context;

    /**
     * the locus dataum constructor
     *
     * @param rodData our reference data
     * @param ref     our reference sequence base position
     * @param context the genome context we're in
     */
    public LocusDatum(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        this.rodData = rodData;
        this.ref = ref;
        this.context = context;
    }

    /**
     * return the Reference order data for this position
     *
     * @return
     */
    public List<ReferenceOrderedDatum> getRodData() {
        return rodData;
    }

    /**
     * return the reference base
     *
     * @return a character representing the reference base
     */
    public char getRef() {
        return ref;
    }

    /**
     * get the locus context at the current position
     *
     * @return
     */
    public LocusContext getContext() {
        return context;
    }

    /**
     * gets the current postion in the sequence, which comes
     * free from underlying data types
     *
     * @return our current GenomeLocation
     */
    public GenomeLoc getSequenceLocation() {
        return this.context.getLocation();
    }
}
