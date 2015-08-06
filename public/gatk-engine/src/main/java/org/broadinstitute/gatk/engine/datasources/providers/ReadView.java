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

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;

import java.util.Arrays;
import java.util.Collection;
/**
 * User: hanna
 * Date: May 22, 2009
 * Time: 12:06:54 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A view into the reads that a provider can provide. 
 */
public class ReadView implements View, Iterable<SAMRecord> {
    /**
     * The iterator into the reads supplied by this provider.
     */
    private GATKSAMIterator reads;

    /**
     * Create a new view of the reads given the current data set.
     * @param provider Source for the data.
     */
    public ReadView( ReadShardDataProvider provider ) {
        reads = provider.getReadIterator();
    }

    /**
     * Other reads and loci conflict with this view.
     * @return Array of reads and loci.
     */
    public Collection<Class<? extends View>> getConflictingViews() {
        return Arrays.<Class<? extends View>>asList(ReadView.class, LocusView.class);
    }

    /**
     * Close the view over these reads.  Note that this method closes just
     * the view into the reads, not the reads themselves.
     */
    public void close() {
        // Don't close the reads.  The provider is responsible for this.
        // Just dispose of the pointer.
        reads = null;
    }

    /**
     * Gets an iterator into the reads supplied by this provider.
     * @return Iterator into the reads that this provider covers.
     */
    public GATKSAMIterator iterator() {
        return reads;    
    }
}
