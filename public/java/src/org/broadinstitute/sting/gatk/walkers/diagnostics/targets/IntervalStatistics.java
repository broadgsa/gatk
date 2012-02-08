package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Short one line description of the walker.
 *
 * @author Mauricio Carneiro
 * @since 2/1/12
 */
class IntervalStatistics {
    private final GenomeLoc interval;
    private final ArrayList<IntervalStatisticLocus> loci;

    private final int minimumCoverageThreshold;
    private final int maximumCoverageThreshold;
    private final int minimumMappingQuality;
    private final int minimumBaseQuality;

    private int preComputedTotalCoverage = -1;         // avoids re-calculating the total sum (-1 means we haven't pre-computed it yet)

    private IntervalStatistics(GenomeLoc interval, ArrayList<IntervalStatisticLocus> loci, int minimumCoverageThreshold, int maximumCoverageThreshold, int minimumMappingQuality, int minimumBaseQuality) {
        this.interval = interval;
        this.loci = loci;
        this.minimumCoverageThreshold = minimumCoverageThreshold;
        this.maximumCoverageThreshold = maximumCoverageThreshold;
        this.minimumMappingQuality = minimumMappingQuality;
        this.minimumBaseQuality = minimumBaseQuality;
    }

    public IntervalStatistics(GenomeLoc interval, int minimumCoverageThreshold, int maximumCoverageThreshold, int minimumMappingQuality, int minimumBaseQuality) {
        this(interval, new ArrayList<IntervalStatisticLocus>(interval.size()), minimumCoverageThreshold, maximumCoverageThreshold, minimumMappingQuality, minimumBaseQuality);

        // Initialize every loci (this way we don't have to worry about non-existent loci in the object
        for (int i = 0; i < interval.size(); i++)
            this.loci.add(i, new IntervalStatisticLocus());

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
     * Calculates the callable status of the entire interval
     *
     * @return the callable status of the entire interval
     */
    public CallableStatus callableStatus() {
        long max = -1;
        CallableStatus maxCallableStatus = null;
        HashMap<CallableStatus, Integer> statusCounts = new HashMap<CallableStatus, Integer>(CallableStatus.values().length);

        // initialize the statusCounts with all callable states
        for (CallableStatus key : CallableStatus.values())
            statusCounts.put(key, 0);

        // calculate the callable status for each locus
        for (int i = 0; i < loci.size(); i++) {
            CallableStatus status = callableStatus(i);
            int count = statusCounts.get(status) + 1;
            statusCounts.put(status, count);

            if (count > max) {
                max = count;
                maxCallableStatus = status;
            }
        }

        return maxCallableStatus;
    }

    public void addLocus(GenomeLoc locus, IntervalStatisticLocus locusData) {
        if (!interval.containsP(locus))
            throw new ReviewedStingException(String.format("Locus %s is not part of the Interval", locus));

        int locusIndex = locus.getStart() - interval.getStart();

        loci.add(locusIndex, locusData);
    }

    /**
     * returns the callable status of this locus without taking the reference base into account.
     *
     * @param locusIndex location in the genome to inquire (only one locus)
     * @return the callable status of a locus
     */
    private CallableStatus callableStatus(int locusIndex) {
        if (loci.get(locusIndex).getCoverage() > maximumCoverageThreshold)
            return CallableStatus.EXCESSIVE_COVERAGE;

        if (loci.get(locusIndex).getCoverage() >= minimumCoverageThreshold)
            return CallableStatus.CALLABLE;

        if (loci.get(locusIndex).getRawCoverage() >= minimumCoverageThreshold)
            return CallableStatus.POOR_QUALITY;

        if (loci.get(locusIndex).getRawCoverage() > 0)
            return CallableStatus.LOW_COVERAGE;

        return CallableStatus.NO_COVERAGE;
    }

    private void calculateTotalCoverage() {
        preComputedTotalCoverage = 0;
        for (IntervalStatisticLocus locus : loci)
            preComputedTotalCoverage += locus.getCoverage();
    }

}
