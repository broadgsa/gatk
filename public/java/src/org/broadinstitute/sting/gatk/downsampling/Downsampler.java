/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.downsampling;

import java.util.Collection;
import java.util.List;

/**
 * The basic downsampler API, with no reads-specific operations
 *
 * @author David Roazen
 */
public interface Downsampler<T> {

    /*
     * Submit one item to the downsampler for consideration . Some downsamplers will be able to determine
     * immediately whether the item survives the downsampling process, while others will need to see
     * more items before making that determination.
     */
    public void submit( T item );

    /*
     * Submit a collection of items to the downsampler for consideration.
     */
    public void submit( Collection<T> items );

    /*
     * Are there items that have survived the downsampling process waiting to be retrieved?
     */
    public boolean hasDownsampledItems();

    /*
     * Return (and remove) all items that have survived downsampling and are waiting to be retrieved.
     */
    public List<T> consumeDownsampledItems();

    /*
     * Are there items stored in this downsampler that it doesn't yet know whether they will
     * ultimately survive the downsampling process?
     */
    public boolean hasPendingItems();

    /*
     * Used to tell the downsampler that no more items will be submitted to it, and that it should
     * finalize any pending items.
     */
    public void signalEndOfInput();

    /*
     * Reset the downsampler to a clean state, devoid of any pending/downsampled items or tracked state
     * information.
     */
    public void clear();
}
