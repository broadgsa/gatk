/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.util;

import java.util.List;
import java.util.Arrays;
import java.io.Closeable;
import java.io.IOException;

/**
 * Utility to close things that implement Closeable
 *
 * @author Kathleen Tibbetts
 */
public class CloserUtil {

    /**
     * Calls close() on <code>obj</code> if it implements Closeable
     *
     * @param obj   The potentially closeable object
     */
    public static void close(Object obj) {
        close(Arrays.asList(obj));
    }

    /**
     * Calls close() on all elements of <code>objs</code> that implement Closeable
     *
     * @param objs   A list of potentially closeable objects
     */
    public static void close(List<Object> objs) {
        for (Object o : objs) {
            if (o instanceof Closeable) {
                try {
                    ((Closeable)o).close();
                }
                catch (IOException ioe) {
                    // Do nothing 
                }
            }
        }
    }
}
