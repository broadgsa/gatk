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

package org.broadinstitute.gatk.utils.collections;

import java.util.List;

/**
 * Represent a permutation of a ordered set or list of elements.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public interface Permutation<E> {

    /**
     * Checks whether this permutation is a partial one of the original list.
     *
     * <p>
     *     A partial permutation is one in that no all original elements take part of.
     * </p>
     *
     * @return {@code true} iff this is a partial permutation.
     */
    public boolean isPartial();

    /**
     * Checks whether this is a trivial permutation where the resulting element list is the same as original.
     *
     * @return {@code true} iff the resulting element list is the same as the original.
     */
    public boolean isNonPermuted();

    /**
     * Given an index on the original list, returns the position of tha element in the resulting list.
     *
     * @param fromIndex the query original element index.
     *
     * @throws IllegalArgumentException if {@code fromIndex} is not a valid index within the original list.
     *
     * @return -1 if that element is not part of the result (partial) permutation, otherwise some number between
     *   0 and {@link #toSize()} - 1.
     */
    public int toIndex(final int fromIndex);

    /**
     * Given an index on the resulting list, it gives you the index of that element on the original list.
     * @param toIndex the query resulting list index.
     *
     * @throws IllegalArgumentException if {@code toIndex} is not a valid index, i.e. in [0,{@link #toSize()}-1).
     *
     * @return a value between 0 and {@link #fromSize()} - 1.
     */
    public int fromIndex(final int toIndex);

    /**
     * Length of the original element list.
     *
     * @return 0 or greater.
     */
    public int fromSize();

    /**
     * Length of the resulting element list.
     *
     * @return 0 or greater.
     */
    public int toSize();

    /**
     * Returns an unmodifiable view to the original element list.
     * @return never {@code null}.
     */
    public List<E> fromList();

    /**
     * Returns an unmodifiable view to the original element list.
     *
     * @return never {@code null}.
     */
    public List<E> toList();
}
