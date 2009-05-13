package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.dataSources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.dataSources.shards.ReadShard;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;

import java.io.File;
import java.util.*;

import edu.mit.broad.picard.filter.FilteringIterator;
import edu.mit.broad.picard.filter.SamRecordFilter;

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
 * @author Mark DePristo
 * @version 0.1
 * @date Apr 29, 2009
 * <p/>
 * Class TraverseDuplicates
 * <p/>
 * This class handles traversing lists of duplicate reads in the new shardable style
 */
public class TraverseDuplicates extends TraversalEngine {
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraverseDuplicates.class);

    private final boolean DEBUG = false;

    /**
     * Creates a new, uninitialized TraversalEngine
     *
     * @param reads SAM/BAM file of reads
     * @param ref   Reference file in FASTA format, assumes a .dict file is also available
     * @param rods  Array of reference ordered data sets
     */
    public TraverseDuplicates(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super(reads, ref, rods);
    }

    private List<SAMRecord> readsAtLoc(final SAMRecord read, PushbackIterator<SAMRecord> iter)
    {
        GenomeLoc site = new GenomeLoc(read);
        ArrayList<SAMRecord> l = new ArrayList<SAMRecord>();

        l.add(read);
        for (SAMRecord read2: iter) {
            GenomeLoc site2 = new GenomeLoc(read2);

            // the next read starts too late
            if ( site2.getStart() != site.getStart() ) {
                //System.out.printf("site = %s, site2 = %s%n", site, site2);
                iter.pushback(read2);
                break;
            } else {
                //System.out.printf("Read is a duplicate: %s%n", read.format());
                l.add(read2);
            }
        }

        return l;
    }

    private Pair<List<SAMRecord>, List<SAMRecord>> splitDuplicates(List<SAMRecord> reads) {
        List<SAMRecord> uniques = new ArrayList<SAMRecord>();
        List<SAMRecord> dups = new ArrayList<SAMRecord>();

        // find the first duplicate
        SAMRecord key = null;
        for ( SAMRecord read : reads ) {
            if ( read.getDuplicateReadFlag() ) {
                // this is our key
                key = read;
                if (DEBUG) logger.debug(String.format("Key %s is a duplicate", read.getReadName()));
                break;
            }
        }

        // At this point, there are two possibilities, we have found at least one dup or not
        // if it's a dup, add it to the dups list, otherwise add it to the uniques list 
        if ( key != null ) {
            final GenomeLoc keyLoc = new GenomeLoc(key);
            final GenomeLoc keyMateLoc = new GenomeLoc(key.getMateReferenceIndex(), key.getMateAlignmentStart(), key.getMateAlignmentStart());

            for ( SAMRecord read : reads ) {
                final GenomeLoc readLoc = new GenomeLoc(read);
                final GenomeLoc readMateLoc = new GenomeLoc(read.getMateReferenceIndex(), read.getMateAlignmentStart(), read.getMateAlignmentStart());
                if (DEBUG) logger.debug(String.format("Examining reads at %s vs. %s at %s / %s vs. %s / %s%n", key.getReadName(), read.getReadName(), keyLoc, keyMateLoc, readLoc, readMateLoc));

                // read and key start at the same place, and either the this read and the key
                // share a mate location or the read is flagged as a duplicate
                if ( readLoc.compareTo(keyLoc) == 0 &&
                        ( readMateLoc.compareTo(keyMateLoc) == 0) ||
                          read.getDuplicateReadFlag() ) {
                    // we are at the same position as the dup and have the same mat pos, it's a dup
                    if (DEBUG) logger.debug(String.format("  => Adding read to dups list: %s%n", read));
                    dups.add(read);
                } else {
                    uniques.add(read);
                }
            }
        } else {
            uniques = reads;
        }

        return new Pair<List<SAMRecord>, List<SAMRecord>>(uniques, dups);
    }

    /**
     * Traverse by reads, given the data and the walker
     * @param sum of type T, the return from the walker
     * @param <M> the generic type
     * @param <T> the return type of the reduce function
     * @return
     */
    public <M, T> T actuallyTraverse(DuplicateWalker<M, T> dupWalker,
                                     Iterator<SAMRecord> readIter,
                                     T sum) {
        // while we still have more reads
        // ok, here's the idea.  We get all the reads that start at the same position in the genome
        // We then split the list of reads into sublists of reads:
        //   -> those with the same mate pair position, for paired reads
        //   -> those flagged as unpaired and duplicated but having the same start and end and
        PushbackIterator<SAMRecord> iter = new PushbackIterator<SAMRecord>(readIter);
        for (SAMRecord read: iter) {
            // get the genome loc from the read
            GenomeLoc site = new GenomeLoc(read);
            List<SAMRecord> reads = readsAtLoc(read, iter);
            Pair<List<SAMRecord>, List<SAMRecord>> split = splitDuplicates(reads);
            List<SAMRecord> uniqueReads =  split.getFirst();
            List<SAMRecord> duplicateReads =  split.getSecond();

            logger.debug(String.format("*** TraverseDuplicates.traverse at %s with %d reads has %d unique and %d duplicate reads",
                    site, reads.size(), uniqueReads.size(), duplicateReads.size()));
            if ( reads.size() != uniqueReads.size() + duplicateReads.size() )
                    throw new RuntimeException(String.format("Bug occurred spliting reads [N=%d] at loc %s into unique [N=%d] and duplicates [N=%d], sizes don't match",
                            reads.size(), uniqueReads.size(), duplicateReads.size()));

            // Jump forward in the reference to this locus location
            LocusContext locus = new LocusContext(site, duplicateReads, Arrays.asList(0));

            // update the number of duplicate sets we've seen
            TraversalStatistics.nRecords++;

            // we still have to fix the locus context provider to take care of this problem with > 1 length contexts
            // LocusContext locus = locusProvider.getLocusContext(site);

            byte[] refBases = new byte[0];

            if ( dupWalker.mapUniqueReadsTooP() ) {
                // Send each unique read to the map function
                for ( SAMRecord unique : uniqueReads ) {
                    List<SAMRecord> l = Arrays.asList(unique);
                    sum = mapOne(dupWalker, l, site, refBases, locus, sum);
                }
            }

            if ( duplicateReads.size() > 0 )
                sum = mapOne(dupWalker, duplicateReads, site, refBases, locus, sum);

            printProgress("dups", site);

            if (this.maxReads > 0 && TraversalStatistics.nRecords > this.maxReads) {
                logger.warn(String.format(("Maximum number of duplicate sets encountered, terminating traversal " + TraversalStatistics.nRecords)));
                break;
            }
        }

        return sum;
    }

    /**
     * Class to filter out un-handle-able reads from the stream.  We currently are skipping
     * unmapped reads, non-primary reads, unaligned reads, and duplicate reads.
     */
    public static class duplicateStreamFilterFunc implements SamRecordFilter {
        SAMRecord lastRead = null;
        public boolean filterOut(SAMRecord rec) {
            boolean result = false;
            String why = "";
            if (rec.getReadUnmappedFlag()) {
                TraversalStatistics.nUnmappedReads++;
                result = true;
                why = "Unmapped";
            } else if (rec.getNotPrimaryAlignmentFlag()) {
                TraversalStatistics.nNotPrimary++;
                result = true;
                why = "Not Primary";
            } else if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                TraversalStatistics.nBadAlignments++;
                result = true;
                why = "No alignment start";
            }
            else {
                result = false;
            }

            if (result) {
                TraversalStatistics.nSkippedReads++;
                //System.out.printf("  [filter] %s => %b %s", rec.getReadName(), result, why);
            } else {
                TraversalStatistics.nReads++;
            }
            return result;
        }
    }

    public <M, T> T mapOne(DuplicateWalker<M, T> dupWalker,
                           List<SAMRecord> readSet,
                           GenomeLoc site,
                           byte[] refBases,
                           LocusContext locus,
                           T sum) {
        final boolean keepMeP = dupWalker.filter(site, refBases, locus, readSet);
        if (keepMeP) {
            M x = dupWalker.map(site, refBases, locus, readSet);
            sum = dupWalker.reduce(x, sum);
        }
        return sum;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // old style interface to the system
    //
    // --------------------------------------------------------------------------------------------------------------
    @Override
    public <M,T> T traverse(Walker<M,T> walker, List<GenomeLoc> locations) {
        if ( walker instanceof DuplicateWalker) {
            Walker x = walker;
            DuplicateWalker<?, ?> dupWalker = (DuplicateWalker<?, ?>)x;
            return (T)this.traverseByRead(dupWalker, locations);
        } else {
            throw new IllegalArgumentException("Walker isn't a duplicate walker!");
        }
    }

    /**
     * Should we deleted at the soonist possible opportunity
     */
    public <M, T> Object traverseByRead(DuplicateWalker<M, T> walker, List<GenomeLoc> locations) {
        samReadIter = initializeReads();

        // Initialize the walker
        walker.initialize();

        // Initialize the sum
        FilteringIterator filterIter = new FilteringIterator(samReadIter, new duplicateStreamFilterFunc());
        T sum = actuallyTraverse(walker, filterIter, walker.reduceInit());

        //printOnTraversalDone("reads", sum);
        walker.onTraversalDone(sum);
        return sum;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // new style interface to the system
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * Traverse by reads, given the data and the walker
     * @param walker the walker to execute over
     * @param shard the shard of data to feed the walker
     * @param sum of type T, the return from the walker
     * @param <M> the generic type
     * @param <T> the return type of the reduce function
     * @return
     */
    public <M, T> T traverse(Walker<M, T> walker,
                             Shard shard,
                             ShardDataProvider dataProvider,
                             T sum) {

        logger.debug(String.format("TraverseDuplicates.traverse Genomic interval is %s", ((ReadShard)shard).getSize()));

        if (!(walker instanceof DuplicateWalker))
            throw new IllegalArgumentException("Walker isn't a duplicate walker!");

        DuplicateWalker<M, T> dupWalker = (DuplicateWalker<M, T>) walker;

        // while we still have more reads
        // ok, here's the idea.  We get all the reads that start at the same position in the genome
        // We then split the list of reads into sublists of reads:
        //   -> those with the same mate pair position, for paired reads
        //   -> those flagged as unpaired and duplicated but having the same start and end and

        FilteringIterator filterIter = new FilteringIterator(dataProvider.getReadIterator(), new duplicateStreamFilterFunc());
        PushbackIterator<SAMRecord> iter = new PushbackIterator<SAMRecord>(filterIter);
        return actuallyTraverse(dupWalker, iter, sum);
    }
}