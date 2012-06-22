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

import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.HashSet;
import java.util.Set;

class ThresHolder {
    public static final ThresHolder DEFAULTS = new ThresHolder(20, 20, 5, 700, 20, 50, 0.5, 0.2, 0.5, 0.2, 0.2, 0.5);

    private final int minimumBaseQuality;
    private final int minimumMappingQuality;

    private final int minimumCoverage;
    private final int maximumCoverage;
    private final int minimumMedianDepth;

    private final int maximumInsertSize;

    private final double votePercentageThreshold;
    private final double lowMedianDepthThreshold;
    private final double badMateStatusThreshold;
    private final double coverageStatusThreshold;
    private final double excessiveCoverageThreshold;
    private final double qualityStatusThreshold;

    public ThresHolder(int minimumBaseQuality,
                       int minimumMappingQuality,
                       int minimumCoverage,
                       int maximumCoverage,
                       int minimumMedianDepth,
                       int maximumInsertSize,
                       double votePercentageThreshold,
                       double lowMedianDepthThreshold,
                       double badMateStatusThreshold,
                       double coverageStatusThreshold,
                       double excessiveCoverageThreshold,
                       double qualityStatusThreshold) {
        this.minimumBaseQuality = minimumBaseQuality;
        this.minimumMappingQuality = minimumMappingQuality;
        this.minimumCoverage = minimumCoverage;
        this.maximumCoverage = maximumCoverage;
        this.minimumMedianDepth = minimumMedianDepth;
        this.maximumInsertSize = maximumInsertSize;
        this.votePercentageThreshold = votePercentageThreshold;
        this.lowMedianDepthThreshold = lowMedianDepthThreshold;
        this.badMateStatusThreshold = badMateStatusThreshold;
        this.coverageStatusThreshold = coverageStatusThreshold;
        this.excessiveCoverageThreshold = excessiveCoverageThreshold;
        this.qualityStatusThreshold = qualityStatusThreshold;
    }

    public int getMinimumCoverage() {
        return minimumCoverage;
    }

    public int getMaximumCoverage() {
        return maximumCoverage;
    }

    public int getMinimumMedianDepth() {
        return minimumMedianDepth;
    }

    public int getMaximumInsertSize() {
        return maximumInsertSize;
    }

    public double getVotePercentageThreshold() {
        return votePercentageThreshold;
    }

    public double getLowMedianDepthThreshold() {
        return lowMedianDepthThreshold;
    }

    public double getBadMateStatusThreshold() {
        return badMateStatusThreshold;
    }

    public double getCoverageStatusThreshold() {
        return coverageStatusThreshold;
    }

    public double getExcessiveCoverageThreshold() {
        return excessiveCoverageThreshold;
    }

    public double getQualityStatusThreshold() {
        return qualityStatusThreshold;
    }

    public int getFilteredCoverage(ReadBackedPileup pileup) {
        return pileup.getBaseAndMappingFilteredPileup(minimumBaseQuality, minimumMappingQuality).depthOfCoverage();
    }

    /**
     * Gets the header lines for the VCF writer
     *
     * @return A set of VCF header lines
     */
    public static Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

        // INFO fields for overall data
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        headerLines.add(new VCFInfoHeaderLine("AVG_INTERVAL_DP", 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a loci divided by interval size."));
        headerLines.add(new VCFInfoHeaderLine("Diagnose Targets", 0, VCFHeaderLineType.Flag, "DiagnoseTargets mode"));

        // FORMAT fields for each genotype
        // todo -- find the appropriate VCF constants
        headerLines.add(new VCFFormatHeaderLine("AVG_INTERVAL_DP", 1, VCFHeaderLineType.Float, "Average depth across the interval. Sum of the depth in a loci divided by interval size."));
        headerLines.add(new VCFFormatHeaderLine("Q1", 1, VCFHeaderLineType.Float, "Lower Quartile of depth distribution."));
        headerLines.add(new VCFFormatHeaderLine("MED", 1, VCFHeaderLineType.Float, "Median of depth distribution."));
        headerLines.add(new VCFFormatHeaderLine("Q3", 1, VCFHeaderLineType.Float, "Upper Quartile of depth Distribution."));


        // FILTER fields
        for (CallableStatus stat : CallableStatus.values())
            headerLines.add(new VCFFilterHeaderLine(stat.name(), stat.description));

        return headerLines;
    }

}
