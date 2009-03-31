package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import java.io.Serializable;

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 2:39:05 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public interface SimpleDataSource extends Serializable {

    /**
     * recommend how many data chunks we should be breaking the file into,
     * as a recommendated number.  If not specified (and even if specified)
     * the chunking data source can make decisions to chunk differently.
     *
     * @param chunkCount
     */
    public void chunk(int chunkCount);


}
