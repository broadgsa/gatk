package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.utils.GATKException;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 4:21:58 PM
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
 * @date Apr 6, 2009
 * <p/>
 * Class SimpleDataSourceLoadException
 * <p/>
 * Generate this on a simple data source load exception
 */
public class SimpleDataSourceLoadException extends GATKException {
    public SimpleDataSourceLoadException(String msg) {
        super(msg);
    }

    public SimpleDataSourceLoadException(String msg,Throwable throwable) {
        super(msg,throwable);
    }
}
