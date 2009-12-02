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
    private char[] bases;

    /**
     * Contructor for a simple, windowless reference context.
     * @param locus locus of interest.
     * @param base reference base at that locus.
     */
    public ReferenceContext( GenomeLoc locus, char base ) {
        this( locus, locus, new char[] { base } );
    }

    public ReferenceContext( GenomeLoc locus, GenomeLoc window, char[] bases ) {
        if( !window.containsP(locus) )
            throw new StingException("Invalid locus or window; window does not contain locus");

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
    public char getBase() {
        return bases[(int)(locus.getStart() - window.getStart())];
    }

    /**
     * Get the base at the given locus.
     * @return The base at the given locus from the reference.
     */
    public int getSimpleBase() {
        return BaseUtils.simpleBaseToBaseIndex(getBase());
    }

    /**
     * All the bases in the window currently being examined.
     * @return All bases available.  If the window is of size [0,0], the array will
     *         contain only the base at the given locus.
     */
    public char[] getBases() {
        return bases;
    }
}
