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

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * Short one line description of the walker.
 *
 * @author Mauricio Carneiro
 * @since 2/1/12
 */
class SampleStatistics {
    private final GenomeLoc interval;
    private final ArrayList<LocusStatistics> loci;

    private final int minimumCoverageThreshold;
    private final int maximumCoverageThreshold;
    private final int minimumMappingQuality;
    private final int minimumBaseQuality;

    private int preComputedTotalCoverage = -1;                                                                          // avoids re-calculating the total sum (-1 means we haven't pre-computed it yet)

    private SampleStatistics(GenomeLoc interval, ArrayList<LocusStatistics> loci, int minimumCoverageThreshold, int maximumCoverageThreshold, int minimumMappingQuality, int minimumBaseQuality) {
        this.interval = interval;
        this.loci = loci;
        this.minimumCoverageThreshold = minimumCoverageThreshold;
        this.maximumCoverageThreshold = maximumCoverageThreshold;
        this.minimumMappingQuality = minimumMappingQuality;
        this.minimumBaseQuality = minimumBaseQuality;
    }

    public SampleStatistics(GenomeLoc interval, int minimumCoverageThreshold, int maximumCoverageThreshold, int minimumMappingQuality, int minimumBaseQuality) {
        this(interval, new ArrayList<LocusStatistics>(interval.size()), minimumCoverageThreshold, maximumCoverageThreshold, minimumMappingQuality, minimumBaseQuality);

        // Initialize every loci (this way we don't have to worry about non-existent loci in the object
        for (int i = 0; i < interval.size(); i++)
            this.loci.add(i, new LocusStatistics());

    }

    public long totalCoverage() {
        if (preComputedTotalCoverage < 0)
            calculateTotalCoverage();
        return preComputedTotalCoverage;
    }

    public double averageCoverage() {
        if (preComputedTotalCoverage < 0)
            calculateTotalCoverage();
        return (double) preComputedTotalCoverage / loci.size();
    }

    /**
     * Calculates the callable statuses of the entire interval
     *
     * @return the callable statuses of the entire interval
     */
    public Set<CallableStatus> getCallableStatuses() {

        Map<CallableStatus, Integer> totals = new HashMap<CallableStatus, Integer>(CallableStatus.values().length);

        // initialize map
        for (CallableStatus status : CallableStatus.values())
            totals.put(status, 0);

        // sum up all the callable statuses for each locus
        for (int i = 0; i < interval.size(); i++) {
            for (CallableStatus status : callableStatus(i)) {
                int count = totals.get(status);

                totals.put(status, count + 1);
            }
        }


        Set<CallableStatus> output = new HashSet<CallableStatus>();

        // double to avoid type casting
        double intervalSize = interval.size();

        double coverageStatusThreshold = 0.20;
        if ((totals.get(CallableStatus.NO_COVERAGE) / intervalSize) > coverageStatusThreshold)
            output.add(CallableStatus.NO_COVERAGE);

        if ((totals.get(CallableStatus.LOW_COVERAGE) / intervalSize) > coverageStatusThreshold)
            output.add(CallableStatus.LOW_COVERAGE);

        double excessiveCoverageThreshold = 0.20;
        if ((totals.get(CallableStatus.EXCESSIVE_COVERAGE) / intervalSize) > excessiveCoverageThreshold)
            output.add(CallableStatus.EXCESSIVE_COVERAGE);

        double qualityStatusThreshold = 0.50;
        if ((totals.get(CallableStatus.POOR_QUALITY) / intervalSize) > qualityStatusThreshold)
            output.add(CallableStatus.POOR_QUALITY);

        if (totals.get(CallableStatus.REF_N) > 0)
            output.add(CallableStatus.REF_N);

        if (output.isEmpty()) {
            output.add(CallableStatus.PASS);
        }
        return output;
    }

    /**
     * Adds a locus to the interval wide stats
     *
     * @param locus  The locus given as a GenomeLoc
     * @param pileup The pileup of that locus
     */
    public void addLocus(GenomeLoc locus, ReadBackedPileup pileup) {
        if (!interval.containsP(locus))
            throw new ReviewedStingException(String.format("Locus %s is not part of the Interval", locus));

        // a null pileup means there nothing ot add
        if (pileup != null) {

            int locusIndex = locus.getStart() - interval.getStart();

            int rawCoverage = pileup.depthOfCoverage();
            int coverage = pileup.getBaseAndMappingFilteredPileup(minimumBaseQuality, minimumMappingQuality).depthOfCoverage();

            LocusStatistics locusData = new LocusStatistics(coverage, rawCoverage);

            loci.add(locusIndex, locusData);
        }
    }

    /**
     * returns the callable status of this locus without taking the reference base into account.
     *
     * @param locusIndex location in the genome to inquire (only one locus)
     * @return the callable status of a locus
     */
    private Set<CallableStatus> callableStatus(int locusIndex) {
        LocusStatistics locus = loci.get(locusIndex);

        return locus.callableStatuses(minimumCoverageThreshold, maximumCoverageThreshold);
    }


    private void calculateTotalCoverage() {
        preComputedTotalCoverage = 0;
        for (LocusStatistics locus : loci)
            preComputedTotalCoverage += locus.getCoverage();
    }

}
