package org.broadinstitute.sting.utils.cmdLine;

import java.lang.annotation.*;

/**
 *
 * User: aaron
 * Date: May 8, 2009
 * Time: 7:12:06 PM
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
 * @date May 8, 2009
 * <p/>
 * @interface ArgumentCollection
 * <p/>
 * This object represents an class, that is a collection of arguments.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.FIELD})
public @interface ArgumentCollection {

}
