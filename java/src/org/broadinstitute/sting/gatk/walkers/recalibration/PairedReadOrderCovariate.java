package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Dec 16, 2009
 */

public class PairedReadOrderCovariate implements StandardCovariate{

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize ( final RecalibrationArgumentCollection rac ) { /* do nothing */ }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {
        return read.getReadPairedFlag() ? (read.getFirstOfPairFlag() ? "First_Read" : "Second_Read") : "Not_Paired";
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return str;
    }

}
