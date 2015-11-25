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

package org.broadinstitute.gatk.engine.traversals;

import htsjdk.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.providers.ReadView;
import org.broadinstitute.gatk.utils.iterators.PushbackIterator;
import org.broadinstitute.gatk.engine.walkers.DuplicateWalker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * @author Mark DePristo
 * @version 0.1
 *          <p/>
 *          Class TraverseDuplicates
 *          <p/>
 *          This class handles traversing lists of duplicate reads in the new shardable style
 */
public class TraverseDuplicates<M,T> extends TraversalEngine<M,T,DuplicateWalker<M,T>,ReadShardDataProvider> {
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(TraverseDuplicates.class);

    /** Turn this to true to enable logger.debug output */
    private final boolean DEBUG = false;

    @Override
    public String getTraversalUnits() {
        return "dups";
    }

    private List<GATKSAMRecord> readsAtLoc(final GATKSAMRecord read, PushbackIterator<SAMRecord> iter) {
        GenomeLoc site = engine.getGenomeLocParser().createGenomeLoc(read);
        ArrayList<GATKSAMRecord> l = new ArrayList<GATKSAMRecord>();

        l.add(read);
        for (SAMRecord read2 : iter) {
            GenomeLoc site2 = engine.getGenomeLocParser().createGenomeLoc(read2);

            // the next read starts too late
            if (site2.getStart() != site.getStart()) {
                iter.pushback(read2);
                break;
            } else {
                l.add((GATKSAMRecord) read2);
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
    protected Set<List<GATKSAMRecord>> uniqueReadSets(List<GATKSAMRecord> reads) {
        Set<List<GATKSAMRecord>> readSets = new LinkedHashSet<List<GATKSAMRecord>>();

        // for each read, find duplicates, and either add the read to its duplicate list or start a new one
        for ( GATKSAMRecord read : reads ) {
            List<GATKSAMRecord> readSet = findDuplicateReads(read, readSets);

            if ( readSet == null ) {
                readSets.add(new ArrayList<GATKSAMRecord>(Arrays.asList(read)));    // copy so I can add to the list
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
    protected List<GATKSAMRecord> findDuplicateReads(GATKSAMRecord read, Set<List<GATKSAMRecord>> readSets ) {
        if ( read.getReadPairedFlag() ) {
            // paired
            final GenomeLoc readMateLoc = engine.getGenomeLocParser().createGenomeLoc(read.getMateReferenceName(), read.getMateAlignmentStart(), read.getMateAlignmentStart());

            for (List<GATKSAMRecord> reads : readSets) {
                GATKSAMRecord key = reads.get(0);

                // read and key start at the same place, and either the this read and the key
                // share a mate location or the read is flagged as a duplicate
                if ( read.getAlignmentStart() == key.getAlignmentStart() && key.getReadPairedFlag() && ( key.getDuplicateReadFlag() || read.getDuplicateReadFlag() ) ) {
                    // at least one has to be marked as a duplicate
                    final GenomeLoc keyMateLoc = engine.getGenomeLocParser().createGenomeLoc(key.getMateReferenceName(), key.getMateAlignmentStart(), key.getMateAlignmentStart());
                    if ( readMateLoc.compareTo(keyMateLoc) == 0 ) {
                        // we are at the same position as the dup and have the same mat pos, it's a dup
                        if (DEBUG) logger.debug(String.format("  => Adding read to dups list: %s %d %s vs. %s", read, reads.size(), readMateLoc, keyMateLoc));
                        return reads;
                    }
                }
            }
        } else {
            for (List<GATKSAMRecord> reads : readSets) {
                GATKSAMRecord key = reads.get(0);
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
     *
     * @return the result type T, the product of all the reduce calls
     */
    public T traverse(DuplicateWalker<M, T> walker,
                      ReadShardDataProvider dataProvider,
                      T sum) {
        PushbackIterator<SAMRecord> iter = new PushbackIterator<SAMRecord>(new ReadView(dataProvider).iterator());

        /**
         * while we still have more reads:
         * ok, here's the idea.  We get all the reads that start at the same position in the genome
         * We then split the list of reads into sublists of reads:
         *   -> those with the same mate pair position, for paired reads
         *   -> those flagged as unpaired and duplicated but having the same start and end
         */
        boolean done = walker.isDone();
        for (SAMRecord read : iter) {
            if ( done ) break;
            // get the genome loc from the read
            GenomeLoc site = engine.getGenomeLocParser().createGenomeLoc(read);

            Set<List<GATKSAMRecord>> readSets = uniqueReadSets(readsAtLoc((GATKSAMRecord) read, iter));
            if ( DEBUG ) logger.debug(String.format("*** TraverseDuplicates.traverse at %s with %d read sets", site, readSets.size()));

            // Jump forward in the reference to this locus location
            AlignmentContext locus = new AlignmentContext(site, new ReadBackedPileupImpl(site));

            // update the number of duplicate sets we've seen
            dataProvider.getShard().getReadMetrics().incrementNumIterations();

            // actually call filter and map, accumulating sum
            final boolean keepMeP = walker.filter(site, locus, readSets);
            if (keepMeP) {
                M x = walker.map(site, locus, readSets);
                sum = walker.reduce(x, sum);
            }

            printProgress(site.getStopLocation());
            done = walker.isDone();
        }

        return sum;
    }
}