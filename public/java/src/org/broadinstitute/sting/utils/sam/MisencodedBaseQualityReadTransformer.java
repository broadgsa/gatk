package org.broadinstitute.sting.utils.sam;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Checks for and errors out (or fixes if requested) when it detects reads with base qualities that are not encoded with
 * phred-scaled quality scores.  Q0 == ASCII 33 according to the SAM specification, whereas Illumina encoding starts at
 * Q64.  The idea here is simple: if we are asked to fix the scores then we just subtract 31 from every quality score.
 * Otherwise, we randomly sample reads (for efficiency) and error out if we encounter a qual that's too high.
 */
public class MisencodedBaseQualityReadTransformer extends ReadTransformer {

    private static final int samplingFrequency = 1000;  // sample 1 read for every 1000 encountered
    private static final int encodingFixValue = 31;  // Illumina_64 - PHRED_33

    private boolean disabled;
    private boolean fixQuals;
    private static int currentReadCounter = 0;

    @Override
    public ApplicationTime initializeSub(final GenomeAnalysisEngine engine, final Walker walker) {
        fixQuals = engine.getArguments().FIX_MISENCODED_QUALS;
        disabled = !fixQuals && engine.getArguments().ALLOW_POTENTIALLY_MISENCODED_QUALS;

        return ReadTransformer.ApplicationTime.ON_INPUT;
    }

    @Override
    public boolean enabled() {
        return !disabled;
    }

    @Override
    public GATKSAMRecord apply(final GATKSAMRecord read) {
        if ( fixQuals )
            return fixMisencodedQuals(read);

        checkForMisencodedQuals(read);
        return read;
    }

    protected static GATKSAMRecord fixMisencodedQuals(final GATKSAMRecord read) {
        final byte[] quals = read.getBaseQualities();
        for ( int i = 0; i < quals.length; i++ ) {
            quals[i] -= encodingFixValue;
        }
        read.setBaseQualities(quals);
        return read;
    }

    protected static void checkForMisencodedQuals(final GATKSAMRecord read) {
        // sample reads randomly for checking
        if ( ++currentReadCounter >= samplingFrequency ) {
            currentReadCounter = 0;

            final byte[] quals = read.getBaseQualities();
            for ( final byte qual : quals ) {
                if ( qual > QualityUtils.MAX_REASONABLE_Q_SCORE )
                    throw new UserException.MisencodedBAM(read, "we encountered an extremely high quality score of " + (int)qual);
            }
        }
    }
}
