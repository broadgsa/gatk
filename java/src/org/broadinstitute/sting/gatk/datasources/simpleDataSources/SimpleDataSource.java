package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.datasources.shards.Shard;

import java.io.Serializable;
import java.util.Iterator;


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
 *
 */

/** This class is the interface for all data sources */
public interface SimpleDataSource extends Serializable {


    /**
     * Query the data source for a region of interest, specified by the genome location.
     * The iterator will generate successive calls
     *
     * @param shard the region
     * @return an iterator of the appropriate type, that is limited by the region
     */
    public Iterator seek(Shard shard) throws SimpleDataSourceLoadException;

}
