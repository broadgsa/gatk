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

package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class IntervalStatistics {

    private final Map<String, SampleStatistics> samples;
    private final GenomeLoc interval;

    private int preComputedTotalCoverage = -1;                                                                          // avoids re-calculating the total sum (-1 means we haven't pre-computed it yet)


    public IntervalStatistics(Set<String> samples, GenomeLoc interval, int minimumCoverageThreshold, int maximumCoverageThreshold, int minimumMappingQuality, int minimumBaseQuality) {
        this.interval = interval;
        this.samples = new HashMap<String, SampleStatistics>(samples.size());
        for (String sample : samples)
            this.samples.put(sample, new SampleStatistics(interval, minimumCoverageThreshold, maximumCoverageThreshold, minimumMappingQuality, minimumBaseQuality));
    }

    public SampleStatistics getSample(String sample) {
        return samples.get(sample);
    }

    public void addLocus(AlignmentContext context) {
        ReadBackedPileup pileup = context.getBasePileup();

        for (String sample : samples.keySet())
            getSample(sample).addLocus(context.getLocation(), pileup.getPileupForSample(sample));
    }

    public long totalCoverage() {
        if (preComputedTotalCoverage < 0)
            calculateTotalCoverage();
        return preComputedTotalCoverage;
    }

    public double averageCoverage() {
        if (preComputedTotalCoverage < 0)
            calculateTotalCoverage();
        return (double) preComputedTotalCoverage / interval.size();
    }

    private void calculateTotalCoverage() {
        preComputedTotalCoverage = 0;
        for (SampleStatistics sample : samples.values())
            preComputedTotalCoverage += sample.totalCoverage();
    }

    /**
     * Return the Callable statuses for the interval as a whole
     * todo -- add a voting system for sample flags and add interval specific statuses
     *
     * @return the callable status(es) for the whole interval
     */
    public Set<CallableStatus> callableStatuses() {
        Set<CallableStatus> output = new HashSet<CallableStatus>();

        for (SampleStatistics sample : samples.values())
            output.addAll(sample.getCallableStatuses());

        return output;
    }
}
