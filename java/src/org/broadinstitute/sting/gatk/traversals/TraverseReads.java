package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadView;
import org.broadinstitute.sting.gatk.datasources.providers.ReadReferenceView;
import org.broadinstitute.sting.gatk.datasources.shards.ReadShard;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.IntervalShard;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
public class TraverseReads extends TraversalEngine {
    final ArrayList<String> x = new ArrayList<String>();

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraverseReads.class);


    /**
     * Creates a new, uninitialized TraversalEngine
     *
     * @param reads SAM/BAM file of reads
     * @param ref   Reference file in FASTA format, assumes a .dict file is also available
     * @param rods  Array of reference ordered data sets
     */
    public TraverseReads(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super(reads, ref, rods);
    }


    /**
     * Traverse by reads, given the data and the walker
     * 
     * @param walker the walker to traverse with
     * @param shard the shard, specifying the range of data to iterate over
     * @param dataProvider the provider of the reads data
     * @param sum the value of type T, specified by the walker, to feed to the walkers reduce function
     * @param <M> the map type of the walker
     * @param <T> the reduce type of the walker
     * @return the reduce variable of the read walker
     */
    public <M, T> T traverse(Walker<M, T> walker,
                             Shard shard,
                             ShardDataProvider dataProvider,
                             T sum) {

        if (shard instanceof ReadShard) {
            logger.debug(String.format("TraverseReads.traverse Genomic interval is %s", ((ReadShard) shard).getSize()));
        } else if (shard instanceof IntervalShard) {
            logger.debug(String.format("TraverseReads.traverse Genomic interval is %s", ((IntervalShard) shard).getGenomeLoc()));
        }


        if (!(walker instanceof ReadWalker))
            throw new IllegalArgumentException("Walker isn't a read walker!");

        if( !dataProvider.hasReads() )
            throw new IllegalArgumentException("Unable to traverse reads; no read data is available.");

        ReadWalker<M, T> readWalker = (ReadWalker<M, T>) walker;
        boolean needsReferenceBasesP = WalkerManager.isRequired(walker, DataSource.REFERENCE_BASES);
        
        ReadView reads = new ReadView(dataProvider);
        ReadReferenceView reference = new ReadReferenceView(dataProvider);

        // while we still have more reads
        for (SAMRecord read : reads) {

            // our locus context
            LocusContext locus = null;

            // an array of characters that represent the reference
            char[] refSeq = null;

            if (needsReferenceBasesP && read.getReferenceIndex() >= 0) {
                // get the genome loc from the read
                GenomeLoc site = GenomeLocParser.createGenomeLoc(read);

                // Jump forward in the reference to this locus location
                locus = new LocusContext(site, Arrays.asList(read), Arrays.asList(0));

                // get the array of characters for the reference sequence, since we're a mapped read
                if( dataProvider.hasReference() )
                    refSeq = reference.getReferenceBases( read );
            }

            // update the number of reads we've seen
            TraversalStatistics.nRecords++;

            final boolean keepMeP = readWalker.filter(refSeq, read);
            if (keepMeP) {
                M x = readWalker.map(refSeq, read);
                sum = readWalker.reduce(x, sum);
            }

            printProgress("reads", null);
        }
        return sum;
    }

    /**
     * Temporary override of printOnTraversalDone.
     * TODO: Add some sort of TE.getName() function once all TraversalEngines are ported.
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone( T sum ) {
        printOnTraversalDone( "reads", sum );
    }
}
