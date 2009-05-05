package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.dataSources.shards.ReadShard;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * User: aaron
 * Date: Apr 24, 2009
 * Time: 10:35:22 AM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
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
     * @param walker the walker to execute over
     * @param shard  the shard of data to feed the walker
     * @param sum    of type T, the return from the walker
     * @param <M>    the generic type
     * @param <T>    the return type of the reduce function
     * @return
     */
    public <M, T> T traverse(Walker<M, T> walker,
                             Shard shard,
                             BoundedReadIterator iter,
                             T sum) {

        logger.debug(String.format("TraverseReads.traverse Genomic interval is %s", ((ReadShard) shard).getSize()));

        if (!(walker instanceof ReadWalker))
            throw new IllegalArgumentException("Walker isn't a read walker!");

        ReadWalker<M, T> readWalker = (ReadWalker<M, T>) walker;

        int readCNT = 0;

        // while we still have more reads
        for (SAMRecord read : iter) {

            // our locus context
            LocusContext locus = null;

            if (read.getReferenceIndex() >= 0) {
                // get the genome loc from the read
                GenomeLoc site = new GenomeLoc(read);

                // Jump forward in the reference to this locus location
                locus = new LocusContext(site, Arrays.asList(read), Arrays.asList(0));
            }

            // update the number of reads we've seen
            TraversalStatistics.nRecords++;


            // we still have to fix the locus context provider to take care of this problem with > 1 length contexts
            // LocusContext locus = locusProvider.getLocusContext(site);

            final boolean keepMeP = readWalker.filter(locus, read);
            if (keepMeP) {
                M x = readWalker.map(locus, read);
                sum = readWalker.reduce(x, sum);
            }

            if (locus != null) { printProgress("loci", locus.getLocation()); }
        }
        //System.err.println(TraversalStatistics.nRecords);
        return sum;
    }

}
