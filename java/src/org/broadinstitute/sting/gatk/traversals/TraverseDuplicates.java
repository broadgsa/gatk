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
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.providers.ReadView;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ManagingReferenceOrderedView;
import org.broadinstitute.sting.gatk.datasources.shards.ReadShard;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * @author Mark DePristo
 * @version 0.1
 *          <p/>
 *          Class TraverseDuplicates
 *          <p/>
 *          This class handles traversing lists of duplicate reads in the new shardable style
 */
public class TraverseDuplicates extends TraversalEngine {
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraverseDuplicates.class);

    /** descriptor of the type */
    private static final String DUPS_STRING = "dups";

    /** Turn this to true to enable logger.debug output */
    private final boolean DEBUG = false;

    private List<SAMRecord> readsAtLoc(final SAMRecord read, PushbackIterator<SAMRecord> iter) {
        GenomeLoc site = GenomeLocParser.createGenomeLoc(read);
        ArrayList<SAMRecord> l = new ArrayList<SAMRecord>();

        l.add(read);
        for (SAMRecord read2 : iter) {
            GenomeLoc site2 = GenomeLocParser.createGenomeLoc(read2);

            // the next read starts too late
            if (site2.getStart() != site.getStart()) {
                iter.pushback(read2);
                break;
            } else {
                l.add(read2);
            }
        }

        return l;
    }

    /**
     * Creates a set of lists of reads, where each list contains reads from the same underlying molecule according
     * to their duplicate flag and their (and mate, if applicable) start/end positions.
     *
     * @param reads the list of reads to split into unique molecular samples
     * @return
     */
    protected Set<List<SAMRecord>> uniqueReadSets(List<SAMRecord> reads) {
        Set<List<SAMRecord>> readSets = new HashSet<List<SAMRecord>>();

        // for each read, find duplicates, and either add the read to its duplicate list or start a new one
        for ( SAMRecord read : reads ) {
            List<SAMRecord> readSet = findDuplicateReads(read, readSets);

            if ( readSet == null ) {
                readSets.add(new ArrayList<SAMRecord>(Arrays.asList(read)));    // copy so I can add to the list
            } else {
                readSet.add(read);
            }
        }

        return readSets;
    }

    /**
     * Find duplicate reads for read in the set of unique reads.  This is effective a duplicate marking algorithm,
     * but it relies for safety's sake on the file itself being marked by a true duplicate marking algorithm.  Pair
     * and single-end read aware.
     *
     * @param read
     * @param readSets
     * @return The list of duplicate reads that read is a member of, or null if it's the only one of its kind
     */
    protected List<SAMRecord> findDuplicateReads(SAMRecord read, Set<List<SAMRecord>> readSets ) {
        if ( read.getReadPairedFlag() ) {
            // paired
            final GenomeLoc readMateLoc = GenomeLocParser.createGenomeLoc(read.getMateReferenceIndex(), read.getMateAlignmentStart(), read.getMateAlignmentStart());

            for (List<SAMRecord> reads : readSets) {
                SAMRecord key = reads.get(0);

                // read and key start at the same place, and either the this read and the key
                // share a mate location or the read is flagged as a duplicate
                if ( read.getAlignmentStart() == key.getAlignmentStart() && key.getReadPairedFlag() && ( key.getDuplicateReadFlag() || read.getDuplicateReadFlag() ) ) {
                    // at least one has to be marked as a duplicate
                    final GenomeLoc keyMateLoc = GenomeLocParser.createGenomeLoc(key.getMateReferenceIndex(), key.getMateAlignmentStart(), key.getMateAlignmentStart());
                    if ( readMateLoc.compareTo(keyMateLoc) == 0 ) {
                        // we are at the same position as the dup and have the same mat pos, it's a dup
                        if (DEBUG) logger.debug(String.format("  => Adding read to dups list: %s %d %s vs. %s", read, reads.size(), readMateLoc, keyMateLoc));
                        return reads;
                    }
                }
            }
        } else {
            for (List<SAMRecord> reads : readSets) {
                SAMRecord key = reads.get(0);
                boolean v = (! key.getReadPairedFlag()) && read.getAlignmentStart() == key.getAlignmentStart() && ( key.getDuplicateReadFlag() || read.getDuplicateReadFlag() ) && read.getReadLength() == key.getReadLength();
                //System.out.printf("%s %s %b %b %d %d %d %d => %b%n",
                //        read.getReadPairedFlag(), key.getReadPairedFlag(), read.getDuplicateReadFlag(), key.getDuplicateReadFlag(),
                //        read.getAlignmentStart(), key.getAlignmentStart(), read.getReadLength(), key.getReadLength(), v);
                if ( v ) {
                    //System.out.printf("Returning reads...%n");
                    return reads;
                }
            }
        }

        return null;
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
                //System.out.printf("  [filter] %s => %b", rec.getReadName(), result);
            } else {
                TraversalStatistics.nReads++;
            }
            return result;
        }
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
     * @param sum    of type T, the return from the walker
     * @param <M>    the generic type
     * @param <T>    the return type of the reduce function
     *
     * @return the result type T, the product of all the reduce calls
     */
    public <M, T> T traverse(Walker<M, T> walker,
                             ShardDataProvider dataProvider,
                             T sum) {
        // safety first :-)
        if (!(walker instanceof DuplicateWalker))
            throw new IllegalArgumentException("Walker isn't a duplicate walker!");
        DuplicateWalker<M, T> dupWalker = (DuplicateWalker<M, T>) walker;

        FilteringIterator filterIter = new FilteringIterator(new ReadView(dataProvider).iterator(), new duplicateStreamFilterFunc());
        PushbackIterator<SAMRecord> iter = new PushbackIterator<SAMRecord>(filterIter);

        /**
         * while we still have more reads:
         * ok, here's the idea.  We get all the reads that start at the same position in the genome
         * We then split the list of reads into sublists of reads:
         *   -> those with the same mate pair position, for paired reads
         *   -> those flagged as unpaired and duplicated but having the same start and end
         */
        for (SAMRecord read : iter) {
            // get the genome loc from the read
            GenomeLoc site = GenomeLocParser.createGenomeLoc(read);

            Set<List<SAMRecord>> readSets = uniqueReadSets(readsAtLoc(read, iter));
            if ( DEBUG ) logger.debug(String.format("*** TraverseDuplicates.traverse at %s with %d read sets", site, readSets.size()));

            // Jump forward in the reference to this locus location
            AlignmentContext locus = new AlignmentContext(site, new ReadBackedPileup(site));

            // update the number of duplicate sets we've seen
            TraversalStatistics.nRecords++;

            // actually call filter and map, accumulating sum
            final boolean keepMeP = dupWalker.filter(site, locus, readSets);
            if (keepMeP) {
                M x = dupWalker.map(site, locus, readSets);
                sum = dupWalker.reduce(x, sum);
            }

            printProgress(DUPS_STRING, site);

            if (this.maximumIterations > 0 && TraversalStatistics.nRecords > this.maximumIterations) {
                logger.warn(String.format(("Maximum number of duplicate sets encountered, terminating traversal " + TraversalStatistics.nRecords)));
                break;
            }
        }

        return sum;
    }

    /**
     * Temporary override of printOnTraversalDone.
     *
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone(T sum) {
        printOnTraversalDone(DUPS_STRING, sum);
    }
}