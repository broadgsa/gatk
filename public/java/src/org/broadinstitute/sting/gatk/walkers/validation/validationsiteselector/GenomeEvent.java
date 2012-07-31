/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.validation.validationsiteselector;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.HashMap;
import java.util.List;


public class GenomeEvent implements Comparable {
    final protected GenomeLoc loc;
    /** A set of the alleles segregating in this context */
    final protected List<Allele> alleles;
//    final protected HashMap<String, Object> attributes;

    public GenomeEvent(GenomeLocParser parser, final String contig, final int start, final int stop, final List<Allele> alleles, HashMap<String, Object> attributes) {
        this.loc = parser.createGenomeLoc(contig, start, stop);
        this.alleles = alleles;
//        this.attributes = attributes;
    }

    // Routine to compare two variant contexts (useful to sort collections of vc's).
    // By default, we want to sort first by contig, then by start location

    public GenomeLoc getGenomeLoc() {
        return loc;
    }
    public int compareTo(final Object o) {
        if (!(o instanceof GenomeEvent))
            throw new ReviewedStingException("BUG: comparing variant context with non-VC object");

        GenomeEvent otherEvent = (GenomeEvent)o;

         return loc.compareTo(otherEvent.getGenomeLoc());
    }

    public VariantContext createVariantContextFromEvent() {
        return new VariantContextBuilder("event", loc.getContig(), loc.getStart(), loc.getStop(), alleles)
                .log10PError(0.0).make();

    }
}
