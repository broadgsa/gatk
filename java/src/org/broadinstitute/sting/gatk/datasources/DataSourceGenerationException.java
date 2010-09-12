package org.broadinstitute.sting.gatk.datasources;

import org.broadinstitute.sting.utils.exceptions.GATKException;

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 9:25:49 AM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */

/**
 * This exception is throw when we're unable to generate a data source,
 * most likely due to an incomplete input source list
 */
public class DataSourceGenerationException extends GATKException {
    public DataSourceGenerationException(String message) {
        super(message);
    }
   
    public DataSourceGenerationException(String message, Throwable throwable) {
        super(message, throwable);
    }
}
