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

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.downsampling.Downsampler;
import org.broadinstitute.sting.gatk.downsampling.PassThroughDownsampler;
import org.broadinstitute.sting.gatk.downsampling.ReservoirDownsampler;

import java.util.*;

/**
 * Divides reads by sample and (if requested) does a preliminary downsampling pass with a ReservoirDownsampler.
 *
 * Note: stores reads by sample ID string, not by sample object
 */
class SamplePartitioner {
    private Map<String, Downsampler<SAMRecord>> readsBySample;

    public SamplePartitioner(final LIBSDownsamplingInfo LIBSDownsamplingInfo, final List<String> samples) {
        readsBySample = new HashMap<String, Downsampler<SAMRecord>>(samples.size());
        for ( String sample : samples ) {
            readsBySample.put(sample, createDownsampler(LIBSDownsamplingInfo));
        }
    }

    private Downsampler<SAMRecord> createDownsampler(final LIBSDownsamplingInfo LIBSDownsamplingInfo) {
        return LIBSDownsamplingInfo.isPerformDownsampling()
                ? new ReservoirDownsampler<SAMRecord>(LIBSDownsamplingInfo.getToCoverage())
                : new PassThroughDownsampler<SAMRecord>();
    }

    public void submitRead(SAMRecord read) {
        String sampleName = read.getReadGroup() != null ? read.getReadGroup().getSample() : null;
        if (readsBySample.containsKey(sampleName))
            readsBySample.get(sampleName).submit(read);
    }

    public void doneSubmittingReads() {
        for ( Map.Entry<String, Downsampler<SAMRecord>> perSampleReads : readsBySample.entrySet() ) {
            perSampleReads.getValue().signalEndOfInput();
        }
    }

    public Collection<SAMRecord> getReadsForSample(String sampleName) {
        if ( ! readsBySample.containsKey(sampleName) )
            throw new NoSuchElementException("Sample name not found");

        return readsBySample.get(sampleName).consumeFinalizedItems();
    }

    public void reset() {
        for ( Map.Entry<String, Downsampler<SAMRecord>> perSampleReads : readsBySample.entrySet() ) {
            perSampleReads.getValue().clear();
            perSampleReads.getValue().reset();
        }
    }
}
