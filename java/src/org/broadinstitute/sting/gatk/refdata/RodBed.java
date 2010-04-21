/*
 * Copyright (c) 2009 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * ﬁles (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, sub ject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata;

import java.util.*;
import java.io.IOException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import net.sf.samtools.util.CloseableIterator;


/**
 *  Simple bed format parser:
 *
 * http://genome.ucsc.edu/FAQ/FAQformat.html
 *
 *
 * User: mdepristo
 * Date: April 20, 2010
 * Time: 10:47:14 AM
 */
public class RodBed extends BasicReferenceOrderedDatum {
    protected GenomeLoc loc;
    private List<String> fields;

    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public RodBed(final String name) {
        super(name);
    }


    // ----------------------------------------------------------------------
    //
    // ROD accessors
    //
    // ----------------------------------------------------------------------
    public GenomeLoc getLocation() {
        return loc;
    }
//
//    public ArrayList<String> getHeader() {
//        return header;
//    }

    public List<String> getFields(final Object key) {
        return fields;
    }

    // ----------------------------------------------------------------------
    //
    // formatting
    //
    // ----------------------------------------------------------------------
    public String toString() { return "BED: " + Utils.join("\t", fields); }

    /**
     * Used by ROD management system to set the data in this ROD associated with a line in a rod
     *
     * @param headerObj
     * @param parts
     * @return
     * @throws IOException
     */
    public boolean parseLine(final Object headerObj, final String[] parts) throws IOException {
        if ( parts.length < 4 )
            throw new StingException("BED format requires at least 3 fields: contig start and stop");

        String contig = parts[0];
        int start = Integer.valueOf(parts[1]) + 1; // 1 indel
        int stop = Integer.valueOf(parts[2]) + 1; // 1 indel
        loc = GenomeLocParser.parseGenomeLoc(contig, start, stop);

        fields = Arrays.asList(parts);
//        for ( int i = 0; i < parts.length; i++ ) {
//            fields.add(parts[i]);
//        }

        return true;
    }
}
