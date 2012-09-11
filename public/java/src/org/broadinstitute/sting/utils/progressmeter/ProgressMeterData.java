package org.broadinstitute.sting.utils.progressmeter;

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
