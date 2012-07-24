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

import java.util.HashSet;
import java.util.Set;

class LocusStatistics {
    private final int coverage;
    private final int rawCoverage;

    public LocusStatistics() {
        this.coverage = 0;
        this.rawCoverage = 0;
    }

    public LocusStatistics(int coverage, int rawCoverage) {
        this.coverage = coverage;
        this.rawCoverage = rawCoverage;
    }

    public int getCoverage() {
        return coverage;
    }

    public int getRawCoverage() {
        return rawCoverage;
    }

    /**
     * Generates all applicable statuses from the coverages in this locus
     *
     * @param thresholds the class contains the statistical threshold for making calls
     * @return a set of all statuses that apply
     */
    public Set<CallableStatus> callableStatuses(ThresHolder thresholds) {
        Set<CallableStatus> output = new HashSet<CallableStatus>();

        // if too much coverage
        if (getCoverage() > thresholds.getMaximumCoverage())
            output.add(CallableStatus.EXCESSIVE_COVERAGE);

        // if not enough coverage
        if (getCoverage() < thresholds.getMinimumCoverage()) {
            // was there a lot of low Qual coverage?
            if (getRawCoverage() >= thresholds.getMinimumCoverage())
                output.add(CallableStatus.POOR_QUALITY);
                // no?
            else {
                // is there any coverage?
                if (getRawCoverage() > 0)
                    output.add(CallableStatus.LOW_COVERAGE);
                else
                    output.add(CallableStatus.COVERAGE_GAPS);
            }
        }

        return output;
    }
}
