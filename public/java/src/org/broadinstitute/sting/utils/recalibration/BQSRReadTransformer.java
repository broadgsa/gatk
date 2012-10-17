package org.broadinstitute.sting.utils.recalibration;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * A ReadTransformer that applies BQSR on the fly to reads
 *
 * User: rpoplin
 * Date: 2/13/12
 */
public class BQSRReadTransformer extends ReadTransformer {
    private boolean enabled;
    private BaseRecalibration bqsr;

    @Override
    public ApplicationTime initializeSub(final GenomeAnalysisEngine engine, final Walker walker) {
        this.enabled = engine.hasBaseRecalibration();
        this.bqsr = engine.getBaseRecalibration();
        final BQSRMode mode = WalkerManager.getWalkerAnnotation(walker, BQSRMode.class);
        return mode.ApplicationTime();
    }

    @Override
    public boolean enabled() {
        return enabled;
    }

    /**
     * initialize a new BQSRReadTransformer that applies BQSR on the fly to incoming reads.
     */
    @Override
    public GATKSAMRecord apply(GATKSAMRecord read) {
        bqsr.recalibrateRead(read);
        return read;
    }
}
