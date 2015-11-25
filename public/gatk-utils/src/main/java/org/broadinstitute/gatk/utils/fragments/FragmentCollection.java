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

package org.broadinstitute.gatk.utils.fragments;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Useful helper class to represent the results of the reads -> fragment calculation.
 *
 * Contains singleton -- objects whose underlying reads do not overlap their mate pair
 * Contains overlappingPairs -- objects whose underlying reads do overlap their mate pair
 *
 * User: ebanks, depristo
 * Date: Jan 10, 2011
 */
public class FragmentCollection<T> {
    Collection<T> singletons;
    Collection<List<T>> overlappingPairs;

    public FragmentCollection(final Collection<T> singletons, final Collection<List<T>> overlappingPairs) {
        this.singletons = singletons == null ? Collections.<T>emptyList() : singletons;
        this.overlappingPairs = overlappingPairs == null ? Collections.<List<T>>emptyList() : overlappingPairs;
    }

    /**
     * Gets the T elements not containing overlapping elements, in no particular order
     *
     * @return
     */
    public Collection<T> getSingletonReads() {
        return singletons;
    }

    /**
     * Gets the T elements containing overlapping elements, in no particular order
     *
     * @return
     */
    public Collection<List<T>> getOverlappingPairs() {
        return overlappingPairs;
    }
}
