/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.jna.clibrary;

import com.sun.jna.Platform;

/**
 * Collection of functions that are in the standard CLibrary but are associated with different headers on different platforms.
 */
public class JNAUtils {
    /**
     * Defined in different places on different systems, this is currently 256 on mac and 64 everywhere else.
     */
    public static final int MAXHOSTNAMELEN;

    /**
     * Maximum path length.
     */
    public static final int MAXPATHLEN = 1024;

    static {
      int maxhostnamelen = 64;
      if (Platform.isMac())
         maxhostnamelen = 256;
      MAXHOSTNAMELEN = maxhostnamelen;
    }

    /**
     * Converts a non-zero int to true, otherwise false.
     * @param val int to check.
     * @return true if val is non-zero.
     */
    public static boolean toBoolean(int val) {
        return val != 0;
    }
}
