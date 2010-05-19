/*
 * Copyright (c) 2009 The Broad Institute
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

package org.broadinstitute.sting.gatk.contexts;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;

/**
 * The section of the reference that overlaps with the given
 * read / locus. 
 *
 * @author hanna
 * @version 0.1
 */

public class ReferenceContext {
    /**
     * The locus.
     */
    private GenomeLoc locus;

    /**
     * The window of reference information around the current locus.
     */
    private GenomeLoc window;

    /**
     * The bases in the window around the current locus.
     */
    private byte[] bases;

    private char[] basesAsCharCached = null;

    /**
     * Contructor for a simple, windowless reference context.
     * @param locus locus of interest.
     * @param base reference base at that locus.
     */
    public ReferenceContext( GenomeLoc locus, byte base ) {
        this( locus, locus, new byte[] { base } );
    }

    public ReferenceContext( GenomeLoc locus, GenomeLoc window, byte[] bases ) {
  //      if( !window.containsP(locus) )
  //          throw new StingException("Invalid locus or window; window does not contain locus");

        this.locus = locus;
        this.window = window;
        this.bases = bases;
    }

    /**
     * The locus currently being examined.
     * @return The current locus.
     */
    public GenomeLoc getLocus() {
        return locus;
    }

    public GenomeLoc getWindow() {
        return window;
    }

    /**
     * Get the base at the given locus.
     * @return The base at the given locus from the reference.
     */
    public byte getBase() {
        return bases[(int)(locus.getStart() - window.getStart())];
    }

    @Deprecated
    public char getBaseAsChar() {
        return (char)getBase();
    }

    /**
     * Get the base at the given locus.
     * @return The base at the given locus from the reference.
     */
    public int getBaseIndex() {
        return BaseUtils.simpleBaseToBaseIndex(getBase());
    }

    /**
     * All the bases in the window currently being examined.
     * @return All bases available.  If the window is of size [0,0], the array will
     *         contain only the base at the given locus.
     */
    public byte[] getBases() {
        return bases;
    }

    @Deprecated
    public char[] getBasesAsChars() {
        if ( basesAsCharCached == null )
            basesAsCharCached = new String(bases).toCharArray();
        return basesAsCharCached;
    }


    /** Extracts from the current window and returns n bases starting at this context's locus (NOT
     * from the window start!). The returned array of chars is newly allocated. If n is too large (runs beyond
     * the right boundary of this context's window), an exception will be thrown. If n==(-1), all bases starting
     * from this context's locus through the end of the window will be returned.
     * @param n number of requested bases including and starting from the current locus
     * @return
     */
    public byte[] getBasesAtLocusAsByte(int n) {

        int start = (int)(locus.getStart()-window.getStart());
        int stop = ( n==(-1) ? bases.length : start+n );

        byte[] b = new byte[stop-start];

        if ( stop > bases.length )
            throw new StingException("Bases beyond the current window requested: window="+window+", requested="+n);

        int i = 0;
        for ( int j = start ;  j < stop ; j++) b[i++]=bases[j];
        return b;
    }

    @Deprecated
    public char[] getBasesAtLocus(int n) {
        return new String(getBasesAtLocusAsByte(n)).toCharArray();
    }
}
