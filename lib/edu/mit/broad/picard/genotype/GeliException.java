/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.genotype;

import edu.mit.broad.picard.PicardException;

/**
 * Generic exception thrown by GELI format machinery.
 *
 * @author Doug Voet
 */
public class GeliException extends PicardException {

    public GeliException(String message, Throwable throwable) {
        super(message, throwable);
    }

    public GeliException(String message) {
        super(message);
    }

}
