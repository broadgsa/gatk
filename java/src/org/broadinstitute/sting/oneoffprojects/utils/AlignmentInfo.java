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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.utils;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Aug 3, 2010
 * Time: 4:10:42 PM
 * To change this template use File | Settings | File Templates.
 */

/** A simple utility class that encapsulates information about a single alignment (offset, strand, overlap, mismatch count).
 *
 */
public class AlignmentInfo {
    private int offset;
    private int overlap;
    private int mm ;
    private Assembly a = null;

    private static int RIDICULOUSLY_LARGE_NUMBER = 1000000000;

    public AlignmentInfo() {
        offset = 0;
        mm = RIDICULOUSLY_LARGE_NUMBER;
        overlap = -1;
        a = null;
    }

    public AlignmentInfo(int mm, int offset, boolean isRc, int overlap, Assembly a) {
        this.offset = (isRc ? (-offset-1) : offset );
        this.overlap = overlap;
        this.mm = mm;
        this.a = a;
    }

    boolean isAligned() { return mm < RIDICULOUSLY_LARGE_NUMBER; }

    public boolean isNegativeStrand() { return offset < 0; }
    public Assembly getAssembly() { return a; }
    public int getOffset() { return ( offset < 0 ? (-offset-1) : offset ); }
    public int getMismatchCount() { return mm; }
    public int getOverlap() { return overlap; }
    public double getMismatchRate() { return isAligned() ? ((double)mm)/overlap : 1.0; }
}
