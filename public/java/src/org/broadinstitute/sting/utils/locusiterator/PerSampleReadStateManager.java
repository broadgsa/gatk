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

package org.broadinstitute.sting.utils.locusiterator;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.CigarOperator;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.downsampling.Downsampler;
import org.broadinstitute.sting.gatk.downsampling.LevelingDownsampler;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * ReadStateManager for a single sample
 *
 * User: depristo
 * Date: 1/13/13
 * Time: 12:28 PM
 */
final class PerSampleReadStateManager implements Iterable<AlignmentStateMachine> {
    private final static Logger logger = Logger.getLogger(ReadStateManager.class);
    private final static boolean CAPTURE_DOWNSAMPLING_STATS = true;

    private List<LinkedList<AlignmentStateMachine>> readStatesByAlignmentStart = new LinkedList<LinkedList<AlignmentStateMachine>>();
    private final Downsampler<LinkedList<AlignmentStateMachine>> levelingDownsampler;
    private int thisSampleReadStates = 0;

    private final int downsamplingTarget;
    private int nSitesNeedingDownsampling = 0;
    private int nSites = 0;

    public PerSampleReadStateManager(final LIBSDownsamplingInfo LIBSDownsamplingInfo) {
        this.downsamplingTarget = LIBSDownsamplingInfo.isPerformDownsampling() ? LIBSDownsamplingInfo.getToCoverage() : -1;
        this.levelingDownsampler = LIBSDownsamplingInfo.isPerformDownsampling()
                ? new LevelingDownsampler<LinkedList<AlignmentStateMachine>, AlignmentStateMachine>(LIBSDownsamplingInfo.getToCoverage())
                : null;
    }

    /**
     * Assumes it can just keep the states linked lists without making a copy
     * @param states the new states to add to this manager
     * @return The change in the number of states, after including states and potentially downsampling
     */
    @Requires("states != null")
    @Ensures("result >= 0")
    public int addStatesAtNextAlignmentStart(LinkedList<AlignmentStateMachine> states) {
        if ( states.isEmpty() ) {
            return 0;
        }

        readStatesByAlignmentStart.add(states);
        int nStatesAdded = states.size();

        if ( isDownsampling() ) {
            captureDownsamplingStats();
            levelingDownsampler.submit(readStatesByAlignmentStart);
            levelingDownsampler.signalEndOfInput();

            nStatesAdded -= levelingDownsampler.getNumberOfDiscardedItems();

            // use returned List directly rather than make a copy, for efficiency's sake
            readStatesByAlignmentStart = levelingDownsampler.consumeFinalizedItems();
            levelingDownsampler.reset();
        }

        thisSampleReadStates += nStatesAdded;
        return nStatesAdded;
    }

    private boolean isDownsampling() {
        return levelingDownsampler != null;
    }

    private AlignmentStateMachine getFirst() {
        if (readStatesByAlignmentStart.isEmpty())
            return null;
        else
            return readStatesByAlignmentStart.get(0).getFirst();
    }

    @Requires("isDownsampling()")
    private void captureDownsamplingStats() {
        if ( CAPTURE_DOWNSAMPLING_STATS ) {
            nSites++;
            final int loc = getFirst().getGenomePosition();
            String message = "Pass through";
            final boolean downsampling = thisSampleReadStates > downsamplingTarget;
            if ( downsampling ) {
                nSitesNeedingDownsampling++;
                message = "Downsampling";
            }

            if ( downsampling || nSites % 10000 == 0 )
                logger.info(String.format("%20s at %s: coverage=%d, max=%d, fraction of downsampled sites=%.2e",
                        message, loc, thisSampleReadStates, downsamplingTarget, (1.0 * nSitesNeedingDownsampling / nSites)));
        }
    }

    /**
     * Is there at least one alignment for this sample in this manager?
     * @return true if there's at least one alignment, false otherwise
     */
    public boolean isEmpty() {
        return readStatesByAlignmentStart.isEmpty();
    }

    public AlignmentStateMachine peek() {
        return isEmpty() ? null : readStatesByAlignmentStart.get(0).peek();
    }

    /**
     * Get the number of read states currently in this manager
     * @return the number of read states
     */
    @Ensures("result >= 0")
    public int size() {
        return thisSampleReadStates;
    }

    /**
     * Advances all read states forward by one element, removing states that are
     * no long aligned to the current position.
     * @return the number of states we're removed after advancing
     */
    public int updateReadStates() {
        int nRemoved = 0;
        final Iterator<AlignmentStateMachine> it = iterator();
        while (it.hasNext()) {
            final AlignmentStateMachine state = it.next();
            final CigarOperator op = state.stepForwardOnGenome();
            if (op == null) {
                // we discard the read only when we are past its end AND indel at the end of the read (if any) was
                // already processed. Keeping the read state that returned null upon stepForwardOnGenome() is safe
                // as the next call to stepForwardOnGenome() will return null again AND will clear hadIndel() flag.
                it.remove();                                                // we've stepped off the end of the object
                nRemoved++;
            }
        }

        return nRemoved;
    }

    // todo -- reimplement
    public Iterator<AlignmentStateMachine> iterator() {
        return new Iterator<AlignmentStateMachine>() {
            private final Iterator<LinkedList<AlignmentStateMachine>> alignmentStartIterator = readStatesByAlignmentStart.iterator();
            private LinkedList<AlignmentStateMachine> currentPositionReadStates;
            private Iterator<AlignmentStateMachine> currentPositionReadStatesIterator;

            @Override
            public boolean hasNext() {
                return  alignmentStartIterator.hasNext() ||
                        (currentPositionReadStatesIterator != null && currentPositionReadStatesIterator.hasNext());
            }

            @Override
            public AlignmentStateMachine next() {
                if ( currentPositionReadStatesIterator == null || ! currentPositionReadStatesIterator.hasNext() ) {
                    currentPositionReadStates = alignmentStartIterator.next();
                    currentPositionReadStatesIterator = currentPositionReadStates.iterator();
                }

                return currentPositionReadStatesIterator.next();
            }

            @Override
            public void remove() {
                currentPositionReadStatesIterator.remove();
                thisSampleReadStates--;

                if ( currentPositionReadStates.isEmpty() ) {
                    alignmentStartIterator.remove();
                }
            }
        };
    }
}
