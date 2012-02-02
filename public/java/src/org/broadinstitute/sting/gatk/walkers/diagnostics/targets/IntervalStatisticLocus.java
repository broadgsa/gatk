package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

/**
 * The definition of a locus for the DiagnoseTargets walker statistics calculation
 *
 * @author Mauricio Carneiro
 * @since 2/3/12
 */
class IntervalStatisticLocus {
    private final byte[] mappingQuality;
    private final byte[] baseQuality;
    private final int coverage;
    private final int rawCoverage;

    public IntervalStatisticLocus(byte[] mappingQuality, byte[] baseQuality, int coverage, int rawCoverage) {
        this.mappingQuality = mappingQuality;
        this.baseQuality = baseQuality;
        this.coverage = coverage;
        this.rawCoverage = rawCoverage;
    }

    public IntervalStatisticLocus() {
        this(new byte[1], new byte[1], 0, 0);
    }

    public int getCoverage() {
        return coverage;
    }

    public int getRawCoverage() {
        return rawCoverage;
    }

}
