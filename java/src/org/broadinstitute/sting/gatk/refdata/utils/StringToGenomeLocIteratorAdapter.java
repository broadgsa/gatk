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

package org.broadinstitute.sting.gatk.refdata.utils;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.bed.BedParser;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
* User: asivache
* Date: Jun 11, 2010
* Time: 2:25:42 PM
* To change this template use File | Settings | File Templates.
*/

/**
 * Wrap this adapter around Iterator<String> to get Iterator<GenomLoc>. Each string coming from the underlying
 * iterator is parsed and converted to GenomeLoc on the fly and the latter is returned on each call to next().
 * Two string formats are currently supported: BED and GATK. This iterator will throw an exception if it fails
 * to parse a string.
 */
public class StringToGenomeLocIteratorAdapter implements Iterator<GenomeLoc> {
    private Iterator<String> it = null;

    public enum FORMAT { BED, GATK };

    FORMAT myFormat = FORMAT.GATK;

    public StringToGenomeLocIteratorAdapter(Iterator<String> it, FORMAT format) {
        this.it = it;
        myFormat = format;
    }

    public StringToGenomeLocIteratorAdapter(Iterator<String> it ) {
        this(it,FORMAT.GATK);
    }

    public boolean hasNext() { return it.hasNext(); }

    public GenomeLoc next() {
        if ( myFormat == FORMAT.GATK ) return GenomeLocParser.parseGenomeInterval( it.next() );
        return BedParser.parseLocation( it.next() );
    }

    public void remove() {
        throw new UnsupportedOperationException("method 'remove' is not supported by this iterator");
    }
}
