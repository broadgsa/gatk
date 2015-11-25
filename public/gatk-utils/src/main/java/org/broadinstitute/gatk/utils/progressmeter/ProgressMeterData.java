/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.progressmeter;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

/**
 * a snapshot of our performance, suitable for storage and later analysis
 */
class ProgressMeterData {
    private final double elapsedSeconds;
    private final long unitsProcessed;
    private final long bpProcessed;

    @Requires({"unitsProcessed >= 0", "bpProcessed >= 0", "elapsedSeconds >= 0"})
    public ProgressMeterData(double elapsedSeconds, long unitsProcessed, long bpProcessed) {
        this.elapsedSeconds = elapsedSeconds;
        this.unitsProcessed = unitsProcessed;
        this.bpProcessed = bpProcessed;
    }

    @Ensures("result >= 0.0")
    public double getElapsedSeconds() {
        return elapsedSeconds;
    }

    @Ensures("result >= 0")
    public long getUnitsProcessed() {
        return unitsProcessed;
    }

    @Ensures("result >= 0")
    public long getBpProcessed() {
        return bpProcessed;
    }

    /** How long in seconds to process 1M traversal units? */
    @Ensures("result >= 0.0")
    public double secondsPerMillionElements() {
        return (elapsedSeconds * 1000000.0) / Math.max(unitsProcessed, 1);
    }

    /** How long in seconds to process 1M bp on the genome? */
    @Ensures("result >= 0.0")
    public double secondsPerMillionBP() {
        return (elapsedSeconds * 1000000.0) / Math.max(bpProcessed, 1);
    }

    /** What fraction of the target intervals have we covered? */
    @Requires("targetSize >= 0")
    @Ensures({"result >= 0.0", "result <= 1.0"})
    public double calculateFractionGenomeTargetCompleted(final long targetSize) {
        return (1.0*bpProcessed) / Math.max(targetSize, 1);
    }
}
