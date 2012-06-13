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
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Reservoir Downsampler: Selects n reads out of a stream whose size is not known in advance, with
 * every read in the stream having an equal chance of being selected for inclusion.
 *
 * An implementation of "Algorithm R" from the paper "Random Sampling with a Reservoir" (Jeffrey Scott Vitter, 1985)
 *
 * @author David Roazen
 */
public class ReservoirDownsampler<T extends SAMRecord> implements ReadsDownsampler<T> {

    private ArrayList<T> reservoir;

    private int targetSampleSize;

    private int totalReadsSeen;

    public ReservoirDownsampler ( int targetSampleSize ) {
        if ( targetSampleSize <= 0 ) {
            throw new ReviewedStingException("Cannot do reservoir downsampling with a sample size <= 0");
        }

        this.targetSampleSize = targetSampleSize;
        clear();
    }

    public void submit ( T newRead ) {
        totalReadsSeen++;

        if ( totalReadsSeen <= targetSampleSize ) {
            reservoir.add(newRead);
        }
        else {
            int randomSlot = GenomeAnalysisEngine.getRandomGenerator().nextInt(totalReadsSeen);
            if ( randomSlot < targetSampleSize ) {
                reservoir.set(randomSlot, newRead);
            }
        }
    }

    public void submit ( Collection<T> newReads ) {
        for ( T read : newReads ) {
            submit(read);
        }
    }

    public boolean hasDownsampledItems() {
        return reservoir.size() > 0;
    }

    public List<T> consumeDownsampledItems() {
        List<T> downsampledItems = reservoir;
        clear();
        return downsampledItems;
    }

    public boolean hasPendingItems() {
        return false;
    }

    public void signalEndOfInput() {
        // NO-OP
    }

    public void clear() {
        reservoir = new ArrayList<T>(targetSampleSize);
        totalReadsSeen = 0;
    }

    public boolean requiresCoordinateSortOrder() {
        return false;
    }
}
