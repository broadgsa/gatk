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

package org.broadinstitute.gatk.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.utils.collections.IndexedSet;

import java.util.Collection;

/**
 * Allele list implementation using and indexed-set.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class IndexedAlleleList<A extends Allele> implements AlleleList<A> {

    private final IndexedSet<A> alleles;

    /**
     * Constructs a new empty allele-list
     */
    public IndexedAlleleList() {
        alleles = new IndexedSet<>();
    }

    /**
     * Constructs a new allele-list from an array of alleles.
     *
     * <p>
     *     Repeats in the input array will be ignored (keeping the first one). The order of alleles in the
     *     resulting list is the same as in the natural traversal of the input collection.
     *
     * </p>
     * @param alleles the original allele array
     *
     * @throws java.lang.IllegalArgumentException if {@code alleles} is {@code null} or contains {@code null}s.
     */
    public IndexedAlleleList(final A ... alleles) {
        this.alleles = new IndexedSet<>(alleles);
    }

    /**
     * Constructs a new allele-list from a collection of alleles.
     *
     * <p>
     *     Repeats in the input collection will be ignored (keeping the first one). The order of alleles in the
     *     resulting list is the same as in the natural traversal of the input collection.
     *
     * </p>
     * @param alleles the original allele collection
     *
     * @throws java.lang.IllegalArgumentException if {@code alleles} is {@code null} or contains {@code null}s.
     */
    public IndexedAlleleList(final Collection<A> alleles) {
        this.alleles = new IndexedSet<>(alleles);
    }

    @Override
    public int alleleCount() {
        return alleles.size();
    }

    @Override
    public int alleleIndex(final A allele) {
        return alleles.indexOf(allele);
    }

    @Override
    public A alleleAt(final int index) {
        return alleles.get(index);
    }
}
