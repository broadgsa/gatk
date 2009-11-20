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

package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.refdata.*;

import java.util.*;

/**
 * A window of variants surrounding the current variant being investigated
 *
 * @author ebanks
 * @version 0.1
 */

public class VariantContextWindow {
    /**
     * The variants.
     */
    private LinkedList<Pair<RefMetaDataTracker, RodVCF>> window = new LinkedList<Pair<RefMetaDataTracker, RodVCF>>();
    private int currentContext;

    /**
     * Contructor for a variant context.
     * @param firstVariants  the first set of variants, comprising the right half of the window
     */
    public VariantContextWindow(List<Pair<RefMetaDataTracker, RodVCF>> firstVariants) {
        int windowSize = (firstVariants == null ? 1 : 2 * firstVariants.size() + 1);
        currentContext = (firstVariants == null ? 0 : firstVariants.size());
        window.addAll(firstVariants);
        while ( window.size() < windowSize )
            window.addFirst(null);
    }

    /**
     * The context currently being examined.
     * @return The current context.
     */
    public Pair<RefMetaDataTracker, RodVCF> getContext() {
        return window.get(currentContext);
    }

    /**
     * The maximum number of elements that can be requested on either end of the current context.
     * @return max.
     */
    public int maxWindowElements() {
        return currentContext;
    }

    /**
     * The window around the context currently being examined.
     * @param elementsToLeft number of earlier contexts to return ()
     * @param elementsToRight number of later contexts to return   ()
     * @return The current context window.
     */
    public Pair<RefMetaDataTracker, RodVCF>[] getWindow(int elementsToLeft, int elementsToRight) {
        if ( elementsToLeft > maxWindowElements() || elementsToRight > maxWindowElements() )
            throw new StingException("Too large a window requested");
        if ( elementsToLeft < 0 || elementsToRight < 0 )
            throw new StingException("Window size cannot be negative");        

        Pair[] array = new Pair[elementsToLeft + elementsToRight + 1];
        ListIterator<Pair<RefMetaDataTracker, RodVCF>> iter = window.listIterator(currentContext - elementsToLeft);
        for (int i = 0; i < elementsToLeft + elementsToRight + 1; i++)
            array[i] = iter.next();
        return array;
    }

    /**
     * Move the window along to the next context
     * @param context The new rightmost context
     */
    public void moveWindow(Pair<RefMetaDataTracker, RodVCF> context) {
        window.removeFirst();
        window.addLast(context);
    }
}