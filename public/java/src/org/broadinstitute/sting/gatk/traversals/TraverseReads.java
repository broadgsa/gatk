package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.ReadBasedReferenceOrderedView;
import org.broadinstitute.sting.gatk.datasources.providers.ReadReferenceView;
import org.broadinstitute.sting.gatk.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadView;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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

/**
 * @author aaron
 * @version 1.0
 * @date Apr 24, 2009
 * <p/>
 * Class TraverseReads
 * <p/>
 * This class handles traversing by reads in the new shardable style
 */
public class TraverseReads<M,T> extends TraversalEngine<M,T,ReadWalker<M,T>,ReadShardDataProvider> {    
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraverseReads.class);

    @Override
    protected String getTraversalType() {
        return "reads";
    }

    /**
     * Traverse by reads, given the data and the walker
     *
     * @param walker the walker to traverse with
     * @param dataProvider the provider of the reads data
     * @param sum the value of type T, specified by the walker, to feed to the walkers reduce function
     * @return the reduce variable of the read walker
     */
    public T traverse(ReadWalker<M,T> walker,
                      ReadShardDataProvider dataProvider,
                      T sum) {

        logger.debug(String.format("TraverseReads.traverse Covered dataset is %s", dataProvider));

        if( !dataProvider.hasReads() )
            throw new IllegalArgumentException("Unable to traverse reads; no read data is available.");

        boolean needsReferenceBasesP = WalkerManager.isRequired(walker, DataSource.REFERENCE_BASES);

        ReadView reads = new ReadView(dataProvider);
        ReadReferenceView reference = new ReadReferenceView(dataProvider);

        // get the reference ordered data
        ReadBasedReferenceOrderedView rodView = new ReadBasedReferenceOrderedView(dataProvider);

        boolean done = walker.isDone();
        // while we still have more reads
        for (SAMRecord read : reads) {
            if ( done ) break;
            // ReferenceContext -- the reference bases covered by the read
            ReferenceContext refContext = null;

            // get the array of characters for the reference sequence, since we're a mapped read
            if (needsReferenceBasesP && !read.getReadUnmappedFlag() && dataProvider.hasReference())
                refContext = reference.getReferenceContext(read);

            // update the number of reads we've seen
            ReadMetrics readMetrics = dataProvider.getShard().getReadMetrics();
            readMetrics.incrementNumIterations();

            // if the read is mapped, create a metadata tracker
            ReadMetaDataTracker tracker = (read.getReferenceIndex() >= 0) ? rodView.getReferenceOrderedDataForRead(read) : null;

            final boolean keepMeP = walker.filter(refContext, (GATKSAMRecord) read);
            if (keepMeP) {
                M x = walker.map(refContext, (GATKSAMRecord) read, tracker); // the tracker can be null
                sum = walker.reduce(x, sum);
            }

            GenomeLoc locus = read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX ? null : engine.getGenomeLocParser().createGenomeLoc(read.getReferenceName(),read.getAlignmentStart());
            printProgress(dataProvider.getShard(),locus);
            done = walker.isDone();
        }
        return sum;
    }
}
