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

package org.broadinstitute.sting.gatk.traversals;

import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.datasources.providers.ReadView;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.ReadShard;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 *
 * User: aaron
 * Date: Apr 24, 2009
 * Time: 10:35:22 AM
 */


/**
 * @author Mark DePristo
 * @version 0.1
 * <p/>
 * Class TraverseDuplicates
 * <p/>
 * This class handles traversing lists of duplicate reads in the new shardable style
 */
public class TraverseDuplicates extends TraversalEngine {
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraverseDuplicates.class);

    private final boolean DEBUG = false;

    private List<SAMRecord> readsAtLoc(final SAMRecord read, PushbackIterator<SAMRecord> iter) {
        GenomeLoc site = GenomeLocParser.createGenomeLoc(read);
        ArrayList<SAMRecord> l = new ArrayList<SAMRecord>();

        l.add(read);
        for (SAMRecord read2 : iter) {
            GenomeLoc site2 = GenomeLocParser.createGenomeLoc(read2);

            // the next read starts too late
            if (site2.getStart() != site.getStart()) {
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
        for (SAMRecord read : reads) {
            if (read.getDuplicateReadFlag()) {
                // this is our key
                key = read;
                if (DEBUG) logger.debug(String.format("Key %s is a duplicate", read.getReadName()));
                break;
            }
        }

        // At this point, there are two possibilities, we have found at least one dup or not
        // if it's a dup, add it to the dups list, otherwise add it to the uniques list 
        if (key != null) {
            final GenomeLoc keyLoc = GenomeLocParser.createGenomeLoc(key);
            final GenomeLoc keyMateLoc = GenomeLocParser.createGenomeLoc(key.getMateReferenceIndex(), key.getMateAlignmentStart(), key.getMateAlignmentStart());

            for (SAMRecord read : reads) {
                final GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);
                final GenomeLoc readMateLoc = GenomeLocParser.createGenomeLoc(read.getMateReferenceIndex(), read.getMateAlignmentStart(), read.getMateAlignmentStart());
                if (DEBUG)
                    logger.debug(String.format("Examining reads at %s vs. %s at %s / %s vs. %s / %s%n", key.getReadName(), read.getReadName(), keyLoc, keyMateLoc, readLoc, readMateLoc));

                // read and key start at the same place, and either the this read and the key
                // share a mate location or the read is flagged as a duplicate
                if (readLoc.compareTo(keyLoc) == 0 &&
                        (readMateLoc.compareTo(keyMateLoc) == 0) ||
                        read.getDuplicateReadFlag()) {
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
     *
     * @param sum of type T, the return from the walker
     * @param <M> the generic type
     * @param <T> the return type of the reduce function
     * @param dupWalker our duplicates walker
     * @param readIter our iterator
     *
     * @return the reduce type, T, the final product of all the reduce calls
     */
    private <M, T> T actuallyTraverse(DuplicateWalker<M, T> dupWalker,
                                      Iterator<SAMRecord> readIter,
                                      T sum) {
        /**
         * while we still have more reads:
         * ok, here's the idea.  We get all the reads that start at the same position in the genome
         * We then split the list of reads into sublists of reads:
         *   -> those with the same mate pair position, for paired reads
         *   -> those flagged as unpaired and duplicated but having the same start and end and
         */
        PushbackIterator<SAMRecord> iter = new PushbackIterator<SAMRecord>(readIter);


        for (SAMRecord read : iter) {

            // get the genome loc from the read
            GenomeLoc site = GenomeLocParser.createGenomeLoc(read);

            List<SAMRecord> reads = readsAtLoc(read, iter);
            Pair<List<SAMRecord>, List<SAMRecord>> split = splitDuplicates(reads);
            List<SAMRecord> uniqueReads = split.getFirst();
            List<SAMRecord> duplicateReads = split.getSecond();

            logger.debug(String.format("*** TraverseDuplicates.traverse at %s with %d reads has %d unique and %d duplicate reads",
                    site, reads.size(), uniqueReads.size(), duplicateReads.size()));
            if (reads.size() != uniqueReads.size() + duplicateReads.size())
                throw new RuntimeException(String.format("Bug occurred spliting reads [N=%d] at loc %s into unique [N=%d] and duplicates [N=%d], sizes don't match",
                        reads.size(), site.toString(), uniqueReads.size(), duplicateReads.size()));

            // Jump forward in the reference to this locus location
            LocusContext locus = new LocusContext(site, duplicateReads, Arrays.asList(0));

            // update the number of duplicate sets we've seen
            TraversalStatistics.nRecords++;

            // we still have to fix the locus context provider to take care of this problem with > 1 length contexts
            // LocusContext locus = locusProvider.getLocusContext(site);

            byte[] refBases = new byte[0];

            if (dupWalker.mapUniqueReadsTooP()) {
                // Send each unique read to the map function
                for (SAMRecord unique : uniqueReads) {
                    List<SAMRecord> l = Arrays.asList(unique);
                    sum = mapOne(dupWalker, uniqueReads, l, site, refBases, locus, sum);
                }
            }

            if (duplicateReads.size() > 0)
                sum = mapOne(dupWalker, uniqueReads, duplicateReads, site, refBases, locus, sum);

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
            if (rec.getReadUnmappedFlag()) {
                TraversalStatistics.nUnmappedReads++;
                result = true;
            } else if (rec.getNotPrimaryAlignmentFlag()) {
                TraversalStatistics.nNotPrimary++;
                result = true;
            } else if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                TraversalStatistics.nBadAlignments++;
                result = true;
            } else {
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
                           List<SAMRecord> uniqueReads,
                           List<SAMRecord> duplicateReads,
                           GenomeLoc site,
                           byte[] refBases,
                           LocusContext locus,
                           T sum) {
        final boolean keepMeP = dupWalker.filter(site, refBases, locus, uniqueReads, duplicateReads);
        if (keepMeP) {
            M x = dupWalker.map(site, refBases, locus, uniqueReads, duplicateReads);
            sum = dupWalker.reduce(x, sum);
        }
        return sum;
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // new style interface to the system
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * Traverse by reads, given the data and the walker
     *
     * @param walker the walker to execute over
     * @param shard  the shard of data to feed the walker
     * @param sum    of type T, the return from the walker
     * @param <M>    the generic type
     * @param <T>    the return type of the reduce function
     *
     * @return the result type T, the product of all the reduce calls
     */
    public <M, T> T traverse(Walker<M, T> walker,
                             Shard shard,
                             ShardDataProvider dataProvider,
                             T sum) {

        logger.debug(String.format("TraverseDuplicates.traverse Genomic interval is %s", ((ReadShard) shard).getSize()));

        if (!(walker instanceof DuplicateWalker))
            throw new IllegalArgumentException("Walker isn't a duplicate walker!");

        DuplicateWalker<M, T> dupWalker = (DuplicateWalker<M, T>) walker;

        // while we still have more reads
        // ok, here's the idea.  We get all the reads that start at the same position in the genome
        // We then split the list of reads into sublists of reads:
        //   -> those with the same mate pair position, for paired reads
        //   -> those flagged as unpaired and duplicated but having the same start and end and

        FilteringIterator filterIter = new FilteringIterator(new ReadView(dataProvider).iterator(), new duplicateStreamFilterFunc());
        PushbackIterator<SAMRecord> iter = new PushbackIterator<SAMRecord>(filterIter);
        return actuallyTraverse(dupWalker, iter, sum);
    }


    /**
     * Temporary override of printOnTraversalDone.
     * TODO: Add some sort of TE.getName() function once all TraversalEngines are ported.
     *
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone(T sum) {
        printOnTraversalDone("reads", sum);
    }
}