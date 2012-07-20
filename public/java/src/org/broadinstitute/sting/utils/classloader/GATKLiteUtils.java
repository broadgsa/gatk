/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.classloader;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 *
 * A set of static utility methods for working with the full vs. Lite GATK build
 */
public class GATKLiteUtils {
    /**
     * Constructor access disallowed...static utility methods only!
     */
    private GATKLiteUtils() { }

    /**
     * Utility method to determine whether this is the lite version of the GATK
     */
    public static boolean isGATKLite() {
        if ( isLiteVersion == null ) {
            try {
                Class.forName(DummyProtectedClassName);
                isLiteVersion = false;
            } catch ( ClassNotFoundException e) {
                isLiteVersion = true;
            }
        }
        return isLiteVersion;
    }
    private static final String DummyProtectedClassName = "org.broadinstitute.sting.gatk.DummyProtectedClass";
    private static Boolean isLiteVersion = null;


    /**
     * Utility method to pull out a protected subclass if possible, otherwise it falls back to a public subclass.
     * Important note: the protected classes MUST implement ProtectedPackageSource!
     *
     * @param interfaceClass    the interface class which the target classes implement
     */
    public static Class getProtectedClassIfAvailable(final Class interfaceClass) {
        List<Class<? extends Object>> classes = new PluginManager<Object>(interfaceClass).getPlugins();
        if ( classes.isEmpty() )
            throw new ReviewedStingException("No classes implementing the interface class " + interfaceClass.getSimpleName() + " were found");

        Class result = null;
        for ( Class<? extends Object> c : classes ) {
            if ( ProtectedPackageSource.class.isAssignableFrom(c) ) {
                result = c;
                break;
            }
        }
        if ( result == null )
            result = classes.get(0);

        return result;
    }
}
