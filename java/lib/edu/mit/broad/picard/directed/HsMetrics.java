package edu.mit.broad.picard.directed;

import edu.mit.broad.picard.metrics.MetricBase;

/**
 * The set of metrics captured that are specific to a hybrid selection analysis.
 *
 * @author Tim Fennell
 */
public class HsMetrics extends MetricBase {
    /** The name of the bait set used in the hybrid selection. */
    public String BAIT_SET;

    /** The number of bases in the reference genome used for alignment. */
    public long GENOME_SIZE;

    /** The number of bases which have one or more baits on top of them. */
    public long BAIT_TERRITORY;

    /** The unique number of target bases in the experiment where target is usually exons etc. */
    public long TARGET_TERRITORY;

    /** Target terrirtoy / bait territory.  1 == perfectly efficient, 0.5 = half of baited bases are not target. */
    public double BAIT_DESIGN_EFFICIENCY;

    /** The total number of reads in the SAM or BAM file examine. */
    public int TOTAL_READS;

    /** The number of reads that pass the vendor's filter. */
    public int PF_READS;

    /** The number of PF reads that are not marked as duplicates. */
    public int PF_UNIQUE_READS;

    /** PF reads / total reads.  The percent of reads passing filter. */
    public double PCT_PF_READS;

    /** PF Unique Reads / Total Reads. */
    public double PCT_PF_UQ_READS;

    /** The number of PF reads that are aligned with mapping score > 0 to the reference genome. */
    public int PF_READS_ALIGNED;

    /** PF Reads Aligned / PF Reads. */
    public double PCT_PF_READS_ALIGNED;

    /** The number of bases in the PF aligned reads that are mapped to a reference base. Accounts for clipping and gaps. */
    public int PF_BASES_ALIGNED;

    /** The number of PF aligned bases that mapped to a baited region of the genome. */
    public long ON_BAIT_BASES;

    /** The number of PF aligned bases that mapped to within a fixed interval of a baited region, but not on a baited region. */
    public long NEAR_BAIT_BASES;

    /** The number of PF aligned bases that mapped to neither on or near a bait. */
    public long OFF_BAIT_BASES;

    /** The number of PF aligned bases that mapped to a targetted region of the genome. */
    public long ON_TARGET_BASES;

    /** On+Near Bait Bases / PF Bases Aligned. */
    public double PCT_SELECTED_BASES;

    /** The percentage of aligned PF bases that mapped neither on or near a bait. */
    public double PCT_OFF_BAIT;

    /** The percentage of on+near bait bases that are on as opposed to near. */
    public double ON_BAIT_VS_SELECTED;

    /** The mean coverage of all baits in the experiment. */
    public double MEAN_BAIT_COVERAGE;

    /** The mean coverage of targets that recieved at least coverage depth = 2 at one base. */
    public double MEAN_TARGET_COVERAGE;

    /** The fold by which the baited region has been amplified above genomic background. */
    public double FOLD_ENRICHMENT;

    /** The number of targets that did not reach coverage=2 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /**
     * The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to
     * the mean coverage level in those targets.
     */
    public double FOLD_80_BASE_PENALTY;


    /**
     * Calculates the metrics in this class that can be derived from other metrics in the class.
     */
    public void calculateDerivedMetrics() {
        BAIT_DESIGN_EFFICIENCY = (double) TARGET_TERRITORY / (double) BAIT_TERRITORY;

        PCT_PF_READS         = PF_READS / (double) TOTAL_READS;
        PCT_PF_UQ_READS      = PF_UNIQUE_READS / (double) TOTAL_READS;
        PCT_PF_READS_ALIGNED = PF_READS_ALIGNED / (double) PF_UNIQUE_READS; 

        double denominator   = (ON_BAIT_BASES + NEAR_BAIT_BASES + OFF_BAIT_BASES);

        PCT_SELECTED_BASES   = (ON_BAIT_BASES + NEAR_BAIT_BASES) / denominator;
        PCT_OFF_BAIT         = OFF_BAIT_BASES / denominator;
        ON_BAIT_VS_SELECTED  = ON_BAIT_BASES / (double) (ON_BAIT_BASES + NEAR_BAIT_BASES);
        MEAN_BAIT_COVERAGE   = ON_BAIT_BASES / (double) BAIT_TERRITORY;
        FOLD_ENRICHMENT = (ON_BAIT_BASES/ denominator) / ((double) BAIT_TERRITORY / GENOME_SIZE);
    }
}
