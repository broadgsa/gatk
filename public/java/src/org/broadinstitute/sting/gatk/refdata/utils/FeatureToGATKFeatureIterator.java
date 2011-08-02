/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.utils;

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.Feature;
import org.broad.tribble.iterators.CloseableTribbleIterator;
import org.broadinstitute.sting.utils.GenomeLocParser;


/**
 * 
 * @author aaron 
 * 
 * Class FeatureToGATKFeatureIterator
 *
 * a wrapper on Tribble feature iterators so that they produce GATKFeatures (which produce GenomeLocs)
 */
public class FeatureToGATKFeatureIterator implements CloseableIterator<GATKFeature> {
    private final GenomeLocParser genomeLocParser;
    private final CloseableTribbleIterator<Feature> iterator;
    private final String name;

    public FeatureToGATKFeatureIterator(GenomeLocParser genomeLocParser,CloseableTribbleIterator<Feature> iter, String name) {
        this.genomeLocParser = genomeLocParser;
        this.name = name;
        this.iterator = iter;
    }

    @Override
    public boolean hasNext() {
        return iterator.hasNext();
    }

    @Override
    public GATKFeature next() {
        return new GATKFeature.TribbleGATKFeature(genomeLocParser,iterator.next(),name);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Why does Iterator have this method? We always throw an exception here");
    }

    @Override
    public void close() {
        // we don't close them anymore
    }
}
