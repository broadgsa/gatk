/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
@PartitionBy(PartitionType.READ)
public abstract class ReadWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    public boolean requiresOrderedReads() { return false; }
    
    // Do we actually want to operate on the context?
    /** Must return true for reads that need to be processed. Reads, for which this method return false will
     * be skipped by the engine and never passed to the walker.
     */
    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        // We are keeping all the reads
        return true;
    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public abstract MapType map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker);
}
