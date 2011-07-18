package org.broadinstitute.sting.gatk.datasources.providers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
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

/**
 * User: hanna
 * Date: May 22, 2009
 * Time: 12:36:14 PM
 *
 */

/** Provides access to the reference over a single read. */

public class ReadReferenceView extends ReferenceView {
    /**
     * Create a view of the reference with respect to a single read.
     *
     * @param provider
     */
    public ReadReferenceView( ShardDataProvider provider ) {
        super(provider);
    }

    protected ReferenceContext.ReferenceContextRefProvider getReferenceBasesProvider( GenomeLoc genomeLoc ) {
        return new Provider(genomeLoc);
    }

    public class Provider implements ReferenceContext.ReferenceContextRefProvider {
        GenomeLoc loc;

        public Provider( GenomeLoc loc ) {
            this.loc = loc;
        }

        public byte[] getBases() {
//            System.out.printf("Getting bases for location %s%n", loc);
//            throw new StingException("x");
            return getReferenceBases(loc);
        }
    }

    public ReferenceContext getReferenceContext( SAMRecord read ) {
        GenomeLoc loc = genomeLocParser.createGenomeLoc(read);
//        byte[] bases = super.getReferenceBases(loc);
//        return new ReferenceContext( loc, loc, bases );
        return new ReferenceContext( genomeLocParser, loc, loc, getReferenceBasesProvider(loc) );
    }

}
