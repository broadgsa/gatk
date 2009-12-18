package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 16, 2009
 * Time: 3:22:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class PairedReadOrderCovariate implements ExperimentalCovariate{

    public void initialize (final RecalibrationArgumentCollection rac ) { /* do nothing */ }

    public final Comparable getValue(final SAMRecord read, final int offset) {
        return read.getReadPairedFlag() ? "Not_Paired" : read.getMateUnmappedFlag() ? "Mate_Unmapped" : read.getFirstOfPairFlag() ? "First_Read" : "Second_Read";
    }

    public final Comparable getValue( final String str ) {
        return str.hashCode();
    }

    // Used to estimate the amount space required for the full data HashMap
    public final int estimatedNumberOfBins() {
        return 4;
    }
}
