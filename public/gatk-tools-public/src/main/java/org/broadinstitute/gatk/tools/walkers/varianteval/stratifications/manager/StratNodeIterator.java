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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.manager;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.*;

/**
 * Helper class for creating iterators over all nodes in the stratification tree
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
class StratNodeIterator<T extends Stratifier> implements Iterator<StratNode<T>> {
    Queue<Iterator<StratNode<T>>> iterators = new LinkedList<Iterator<StratNode<T>>>();
    Iterator<StratNode<T>> currentIterator;

    StratNodeIterator(final StratNode<T> root) {
        currentIterator = Collections.singleton(root).iterator();
        for ( final StratNode<T> subNode : root.subnodes.values() )
            iterators.add(new StratNodeIterator<T>(subNode));
    }

    @Override
    public boolean hasNext() {
        return currentIterator.hasNext() || ! iterators.isEmpty();
    }

    @Override
    public StratNode<T> next() {
        if ( currentIterator.hasNext() )
            return currentIterator.next();
        else if ( ! iterators.isEmpty() ) {
            currentIterator = iterators.poll();
            return next();
        } else {
            throw new IllegalStateException("Next called on empty iterator");
        }
    }

    @Override
    public void remove() {
        throw new ReviewedGATKException("Cannot remove from StratNode iterator");
    }
}
