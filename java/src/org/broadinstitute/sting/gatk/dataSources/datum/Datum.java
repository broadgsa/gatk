package org.broadinstitute.sting.gatk.dataSources.datum;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.Serializable;
/**
 *
 * User: aaron
 * Date: Mar 30, 2009
 * Time: 1:32:34 PM
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
 * interface Datum
 * <p/>
 * The interface for all Datum Types.
 */
public interface Datum extends Serializable {

    // this function is used for tracking where we are in a genome
    public GenomeLoc getSequenceLocation();
}
