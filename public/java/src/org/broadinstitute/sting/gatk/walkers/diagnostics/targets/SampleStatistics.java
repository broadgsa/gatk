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
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * The statistics calculator for a specific sample given the interval
 */
class SampleStatistics {
    private final GenomeLoc interval;
    private final ArrayList<LocusStatistics> loci;

    private int[] preSortedDepths = null;
    private int preComputedTotalCoverage = -1;         // avoids re-calculating the total sum (-1 means we haven't pre-computed it yet)

    private int nReads = -1;
    private int nBadMates = -1;

    private SampleStatistics(GenomeLoc interval, ArrayList<LocusStatistics> loci) {
        this.interval = interval;
        this.loci = loci;
        nReads = 0;
        nBadMates = 0;
    }

    public SampleStatistics(GenomeLoc interval) {
        this(interval, new ArrayList<LocusStatistics>(interval.size()));

        // Initialize every loci (this way we don't have to worry about non-existent loci in the object
        for (int i = 0; i < interval.size(); i++)
            this.loci.add(new LocusStatistics());

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
     * Calculates the callable statuses of the entire sample
     *
     * @param thresholds the class contains the statistical threshold for making calls
     * @return the callable statuses of the entire sample
     */
    public Set<CallableStatus> getCallableStatuses(ThresHolder thresholds) {
        // We check if reads are present ot prevent div / 0 exceptions
        if (nReads == 0) {
            return Collections.singleton(CallableStatus.NO_READS);
        }

        Set<CallableStatus> output = new HashSet<CallableStatus>();
        Map<CallableStatus, Double> totals = new HashMap<CallableStatus, Double>(CallableStatus.values().length);

        // initialize map
        for (CallableStatus status : CallableStatus.values())
            totals.put(status, 0.0);

        // sum up all the callable statuses for each locus
        for (int i = 0; i < interval.size(); i++) {
            for (CallableStatus status : callableStatus(i, thresholds)) {
                double count = totals.get(status);

                totals.put(status, count + 1);
            }
        }

        double intervalSize = interval.size();

        if (((double) nBadMates / nReads) >= thresholds.getBadMateStatusThreshold())
            output.add(CallableStatus.BAD_MATE);

        if ((totals.get(CallableStatus.COVERAGE_GAPS) / intervalSize) >= thresholds.getCoverageStatusThreshold())
            output.add(CallableStatus.COVERAGE_GAPS);

        if ((totals.get(CallableStatus.LOW_COVERAGE) / intervalSize) >= thresholds.getCoverageStatusThreshold())
            output.add(CallableStatus.LOW_COVERAGE);

        if ((totals.get(CallableStatus.EXCESSIVE_COVERAGE) / intervalSize) >= thresholds.getExcessiveCoverageThreshold())
            output.add(CallableStatus.EXCESSIVE_COVERAGE);

        if ((totals.get(CallableStatus.POOR_QUALITY) / intervalSize) >= thresholds.getQualityStatusThreshold())
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
     * @param locus      The locus given as a GenomeLoc
     * @param pileup     The pileup of that locus, this exclusively contains the sample
     * @param thresholds the class contains the statistical threshold for making calls
     */
    public void addLocus(GenomeLoc locus, ReadBackedPileup pileup, ThresHolder thresholds) {
        if (!interval.containsP(locus))
            throw new ReviewedStingException(String.format("Locus %s is not part of the Interval %s", locus, interval));

        // a null pileup means there nothing ot add
        if (pileup != null) {

            int locusIndex = locus.getStart() - interval.getStart();

            int rawCoverage = pileup.depthOfCoverage();
            int coverage = thresholds.getFilteredCoverage(pileup);

            LocusStatistics locusData = new LocusStatistics(coverage, rawCoverage);

            loci.set(locusIndex, locusData);

            for (GATKSAMRecord read : pileup.getReads())
                processRead(read, thresholds);
        }
    }

    private void processRead(GATKSAMRecord read, ThresHolder thresholds) {
        // Was this read already processed?
        if (read.getTemporaryAttribute("checkedBadMate") == null) {
            nReads++;
            if (!hasValidMate(read, thresholds))
                nBadMates++;
            read.setTemporaryAttribute("checkedBadMate", true);
        }
    }

    /**
     * returns the callable status of a given locus without taking the reference base into account.
     *
     * @param locusIndex location in the genome to inquire (only one locus)
     * @param thresholds the class contains the statistical threshold for making calls
     * @return the callable status of a locus
     */
    private Set<CallableStatus> callableStatus(int locusIndex, ThresHolder thresholds) {
        LocusStatistics locus = loci.get(locusIndex);

        return locus.callableStatuses(thresholds);
    }

    private void calculateTotalCoverage() {
        preComputedTotalCoverage = 0;
        for (LocusStatistics locus : loci)
            preComputedTotalCoverage += locus.getCoverage();
    }

    public double getQuantileDepth(double percentage) {
        if (preSortedDepths == null)
            getDepthsAsSortedArray();

        return getQuartile(preSortedDepths, percentage);
    }

    static double getQuartile(int[] data, double percentage) {
        int size = data.length;
        if (size == 1)
            return (double) data[0];

        if (percentage == 0.5) {
            return getMedian(data);
        }

        double position = (size - 1.0) / 2;
        if (percentage == 0.25) {
            // if the position is a whole number
            return getMedian(Arrays.copyOfRange(data, 0, (int) position + 1));

        }
        if (percentage == 0.75) {
            if (position % 1 == 0) {
                return getMedian(Arrays.copyOfRange(data, (int) position, size));
            } else {
                return getMedian(Arrays.copyOfRange(data, (int) position + 1, size));
            }
        }
        return -1;
    }

    // Assumes data is sorted
    private static double getMedian(int[] data) {
        double size = (double) data.length;
        if (size == 1)
            return (double) data[0];

        double position = (size - 1.0) / 2;

        if (position % 1 == 0)
            return (double) data[(int) position];

        else {
            double high = (double) data[(int) Math.ceil(position)];
            double low = (double) data[(int) Math.floor(position)];

            return (high + low) / 2;

        }

    }

    private void getDepthsAsSortedArray() {
        preSortedDepths = new int[loci.size()];

        for (int i = 0; i < loci.size(); i++)
            preSortedDepths[i] = loci.get(i).getCoverage();

        Arrays.sort(preSortedDepths);
    }

    boolean hasValidMate(GATKSAMRecord read, ThresHolder thresholds) {
        /** Check the following
         * Does it have a pair?
         * reasonable insert size?
         * inverted?
         * same orientation?
         * same contig?
         * is pair mapped?
         * todo - is forced mate?
         *
         */

        // has NO pair
        if (!read.getReadPairedFlag())
            return false;

        // different contigs
        if (read.getMateReferenceIndex() != read.getReferenceIndex())
            return false;

        // unmapped
        if (read.getMateUnmappedFlag() || read.getReadUnmappedFlag())
            return false;

        // same orientation
        if (read.getReadNegativeStrandFlag() == read.getMateNegativeStrandFlag())
            return false;

        // inverted
        if (read.getReadNegativeStrandFlag() ==
                read.getAlignmentStart() < read.getMateAlignmentStart())
            return false;

        // TODO note: IGV uses a different algorithm for insert size, there should be a common util class that does this for you
        // mates are too far apart
        if (Math.abs(read.getAlignmentStart() - read.getMateAlignmentStart()) > thresholds.getMaximumInsertSize())
            return false;

        return true;
    }

    public int getnReads() {
        return nReads;
    }

    public int getnBadMates() {
        return nBadMates;
    }
}
