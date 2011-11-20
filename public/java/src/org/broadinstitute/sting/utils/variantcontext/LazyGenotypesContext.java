/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.utils.variantcontext;

import org.broadinstitute.sting.utils.codecs.vcf.VCFParser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * [Short one sentence description of this walker]
 * <p/>
 * <p>
 * [Functionality of this walker]
 * </p>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 *  </pre>
 *
 * @author Your Name
 * @since Date created
 */
public class LazyGenotypesContext extends GenotypesContext {
    final VCFParser parser;
    String unparsedGenotypeData;
    final List<Allele> alleles;
    final String contig;
    final int start;
    final int nUnparsedGenotypes;

    boolean loaded = false;

    private final static ArrayList<Genotype> EMPTY = new ArrayList<Genotype>(0);

    public LazyGenotypesContext(final VCFParser parser, final String unparsedGenotypeData,
                                final String contig, final int start, final List<Allele> alleles,
                                int nUnparsedGenotypes ) {
        super(EMPTY, false);
        this.unparsedGenotypeData = unparsedGenotypeData;
        this.start = start;
        this.parser = parser;
        this.contig = contig;
        this.alleles = alleles;
        this.nUnparsedGenotypes = nUnparsedGenotypes;
    }

    @Override
    protected ArrayList<Genotype> getGenotypes() {
        if ( ! loaded ) {
            //System.out.printf("Loading genotypes... %s:%d%n", contig, start);
            GenotypesContext subcontext = parser.createGenotypeMap(unparsedGenotypeData, alleles, contig, start);
            notToBeDirectlyAccessedGenotypes = subcontext.notToBeDirectlyAccessedGenotypes;
            sampleNamesInOrder = subcontext.sampleNamesInOrder;
            sampleNameToOffset = subcontext.sampleNameToOffset;
            cacheIsInvalid = false;
            loaded = true;
            unparsedGenotypeData = null;

            // warning -- this path allows us to create a VariantContext that doesn't run validateGenotypes()
            // That said, it's not such an important routine -- it's just checking that the genotypes
            // are well formed w.r.t. the alleles list, but this will be enforced within the VCFCodec
        }

        return notToBeDirectlyAccessedGenotypes;
    }

    protected synchronized void buildCache() {
        if ( cacheIsInvalid ) {
            getGenotypes(); // will load up all of the necessary data
        }
    }

    @Override
    public boolean isEmpty() {
        // optimization -- we know the number of samples in the unparsed data, so use it here to
        // avoid parsing just to know if the genotypes context is empty
        return loaded ? super.isEmpty() : nUnparsedGenotypes == 0;
    }

    @Override
    public int size() {
        // optimization -- we know the number of samples in the unparsed data, so use it here to
        // avoid parsing just to know the size of the context
        return loaded ? super.size() : nUnparsedGenotypes;
    }

    public String getUnparsedGenotypeData() {
        return unparsedGenotypeData;
    }
}
