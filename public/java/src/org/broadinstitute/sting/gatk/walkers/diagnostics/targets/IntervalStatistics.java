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
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

class IntervalStatistics {

    private final Map<String, SampleStatistics> samples;
    private final GenomeLoc interval;
    private boolean hasNref = false;

    private int preComputedTotalCoverage = -1;                                                                          // avoids re-calculating the total sum (-1 means we haven't pre-computed it yet)

    /*
    private double minMedianDepth = 20.0;
    private double badMedianDepthPercentage = 0.20;
    private double votePercentage = 0.50;
    */
    public IntervalStatistics(Set<String> samples, GenomeLoc interval/*, int minimumCoverageThreshold, int maximumCoverageThreshold, int minimumMappingQuality, int minimumBaseQuality*/) {
        this.interval = interval;
        this.samples = new HashMap<String, SampleStatistics>(samples.size());
        for (String sample : samples)
            this.samples.put(sample, new SampleStatistics(interval /*, minimumCoverageThreshold, maximumCoverageThreshold, minimumMappingQuality, minimumBaseQuality*/));
    }

    public SampleStatistics getSample(String sample) {
        return samples.get(sample);
    }

    public GenomeLoc getInterval() {
        return interval;
    }

    /**
     * The function to populate data into the Statistics from the walker.
     * This takes the input and manages passing the data to the SampleStatistics and Locus Statistics
     *
     * @param context    The alignment context given from the walker
     * @param ref        the reference context given from the walker
     * @param thresholds the class contains the statistical threshold for making calls
     */
    public void addLocus(AlignmentContext context, ReferenceContext ref, ThresHolder thresholds) {
        ReadBackedPileup pileup = context.getBasePileup();

        //System.out.println(ref.getLocus().toString());

        Map<String, ReadBackedPileup> samplePileups = pileup.getPileupsForSamples(samples.keySet());

        for (Map.Entry<String, ReadBackedPileup> entry : samplePileups.entrySet()) {
            String sample = entry.getKey();
            ReadBackedPileup samplePileup = entry.getValue();
            SampleStatistics sampleStatistics = samples.get(sample);

            if (sampleStatistics == null)
                throw new ReviewedStingException(String.format("Trying to add locus statistics to a sample (%s) that doesn't exist in the Interval.", sample));

            sampleStatistics.addLocus(context.getLocation(), samplePileup, thresholds);
        }

        if (!hasNref && ref.getBase() == 'N')
            hasNref = true;
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
     * todo -- add missingness filter
     *
     * @param thresholds the class contains the statistical threshold for making calls
     * @return the callable status(es) for the whole interval
     */
    public Set<CallableStatus> callableStatuses(ThresHolder thresholds) {
        Set<CallableStatus> output = new HashSet<CallableStatus>();

        // Initialize the Map
        Map<CallableStatus, Integer> votes = new HashMap<CallableStatus, Integer>();
        for (CallableStatus status : CallableStatus.values())
            votes.put(status, 0);

        // tally up the votes
        for (SampleStatistics sample : samples.values())
            for (CallableStatus status : sample.getCallableStatuses(thresholds))
                votes.put(status, votes.get(status) + 1);

        // output tall values above the threshold
        for (CallableStatus status : votes.keySet()) {
            if (votes.get(status) > (samples.size() * thresholds.getVotePercentageThreshold()) && !(status.equals(CallableStatus.PASS)))
                output.add(status);
        }


        if (hasNref)
            output.add(CallableStatus.REF_N);

        // get median DP of each sample
        int nLowMedianDepth = 0;
        for (SampleStatistics sample : samples.values()) {
            if (sample.getQuantileDepth(0.5) < thresholds.getMinimumMedianDepth())
                nLowMedianDepth++;
        }

        if (nLowMedianDepth > (samples.size() * thresholds.getLowMedianDepthThreshold()))
            output.add(CallableStatus.LOW_MEDIAN_DEPTH);

        return output;
    }
}
