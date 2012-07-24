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

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Positional Downsampler: When eliminating reads, try to do so evenly based on the alignment start positions
 *
 * @author David Roazen
 */
public class PositionalDownsampler<T extends SAMRecord> implements ReadsDownsampler<T> {

    private int targetCoverage;

    private ReservoirDownsampler<T> reservoir;

    private int currentContigIndex;

    private int currentAlignmentStart;

    private LinkedList<PositionalReadGrouping> pendingReads;

    private ArrayList<T> finalizedReads;

    public PositionalDownsampler ( int targetCoverage ) {
        this.targetCoverage = targetCoverage;
        clear();
    }

    public void submit ( T newRead ) {
        if ( readIsPastCurrentPosition(newRead) ) {
            updateAndDownsamplePendingReads();
        }

        reservoir.submit(newRead);
        updateCurrentPosition(newRead);
    }

    public void submit ( Collection<T> newReads ) {
        for ( T read : newReads ) {
            submit(read);
        }
    }

    public boolean hasDownsampledItems() {
        return finalizedReads.size() > 0;
    }

    public List<T> consumeDownsampledItems() {
        List<T> toReturn = finalizedReads;
        finalizedReads = new ArrayList<T>();
        return toReturn;
    }

    public boolean hasPendingItems() {
        return pendingReads.size() > 0;
    }

    public void signalEndOfInput() {
        updateAndDownsamplePendingReads();

        for ( PositionalReadGrouping group : pendingReads ) {
            group.finalizeAllActiveReads();
            finalizedReads.addAll(group.getFinalizedReads());
        }

        pendingReads.clear();
    }

    public void clear() {
        reservoir = new ReservoirDownsampler<T>(targetCoverage);
        pendingReads = new LinkedList<PositionalReadGrouping>();
        finalizedReads = new ArrayList<T>();
    }

    public boolean requiresCoordinateSortOrder() {
        return true;
    }

    private void updateCurrentPosition ( T read ) {
        currentContigIndex = read.getReferenceIndex();
        currentAlignmentStart = read.getAlignmentStart();
    }

    private boolean readIsPastCurrentPosition ( T read ) {
        return read.getReferenceIndex() != currentContigIndex || read.getAlignmentStart() > currentAlignmentStart;
    }

    private void updateAndDownsamplePendingReads() {
        finalizeOutOfScopeReads();

        List<T> oldLocusReads = reservoir.consumeDownsampledItems();
        pendingReads.add(new PositionalReadGrouping(oldLocusReads, currentContigIndex, currentAlignmentStart));

        downsampleOverlappingGroups();
    }

    private void finalizeOutOfScopeReads() {
        Iterator<PositionalReadGrouping> iter = pendingReads.iterator();
        boolean noPrecedingUnfinalizedGroups = true;

        while ( iter.hasNext() ) {
            PositionalReadGrouping currentGroup = iter.next();
            currentGroup.finalizeActiveReadsBeforePosition(currentContigIndex, currentAlignmentStart);

            if ( currentGroup.isFinalized() && noPrecedingUnfinalizedGroups ) {
                iter.remove();
                finalizedReads.addAll(currentGroup.getFinalizedReads());
            }
            else {
                noPrecedingUnfinalizedGroups = false;
            }
        }
    }

    private void downsampleOverlappingGroups() {
        int[] groupReadCounts = new int[pendingReads.size()];
        int totalCoverage = 0;
        int numActiveGroups = 0;
        int currentGroup = 0;

        for ( PositionalReadGrouping group : pendingReads ) {
            groupReadCounts[currentGroup] = group.numActiveReads();
            totalCoverage += groupReadCounts[currentGroup];

            if ( groupReadCounts[currentGroup] > 0 ) {
                numActiveGroups++;
            }

            currentGroup++;
        }

        if ( totalCoverage <= targetCoverage ) {
            return;
        }

        int numReadsToRemove = Math.min(totalCoverage - targetCoverage, totalCoverage - numActiveGroups);
        currentGroup = 0;

        while ( numReadsToRemove > 0  ) {
            if ( groupReadCounts[currentGroup] > 1 ) {
                groupReadCounts[currentGroup]--;
                numReadsToRemove--;
            }

            currentGroup = (currentGroup + 1) % groupReadCounts.length;
        }

        currentGroup = 0;
        for ( PositionalReadGrouping group : pendingReads ) {
            if ( ! group.isFinalized() ) {
                group.downsampleActiveReads(groupReadCounts[currentGroup]);
            }
            currentGroup++;
        }
    }

    private class PositionalReadGrouping {
        private List<T> activeReads;
        private List<T> finalizedReads;

        private int contig;
        private int alignmentStart;

        public PositionalReadGrouping( Collection<T> reads, int contig, int alignmentStart ) {
            activeReads = new LinkedList<T>(reads);
            finalizedReads = new ArrayList<T>();
            this.contig = contig;
            this.alignmentStart = alignmentStart;
        }

        public int numActiveReads() {
            return activeReads.size();
        }

        public boolean isFinalized() {
            return activeReads.size() == 0;
        }

        public List<T> getFinalizedReads() {
            return finalizedReads;
        }

        public void finalizeActiveReadsBeforePosition( int contig, int position ) {
            if ( this.contig != contig ) {
                finalizeAllActiveReads();
                return;
            }

            Iterator<T> iter = activeReads.iterator();

            while ( iter.hasNext() ) {
                T read = iter.next();
                if ( read.getAlignmentEnd() < position ) {
                    iter.remove();
                    finalizedReads.add(read);
                }
            }
        }

        public void finalizeAllActiveReads() {
            finalizedReads.addAll(activeReads);
            activeReads.clear();
        }

        public void downsampleActiveReads( int numReadsToKeep ) {
            if ( numReadsToKeep > activeReads.size() || numReadsToKeep < 0 ) {
                throw new ReviewedStingException(String.format("Cannot retain %d reads out of %d total reads",
                                                               numReadsToKeep, activeReads.size()));
            }

            BitSet itemsToKeep = new BitSet(activeReads.size());
            for ( Integer selectedIndex : MathUtils.sampleIndicesWithoutReplacement(activeReads.size(), numReadsToKeep) ) {
                itemsToKeep.set(selectedIndex);
            }

            int currentIndex = 0;
            Iterator<T> iter = activeReads.iterator();

            while ( iter.hasNext() ) {
                T read = iter.next();

                if ( ! itemsToKeep.get(currentIndex) ) {
                    iter.remove();
                }

                currentIndex++;
            }
        }

    }
}
