package org.broadinstitute.sting.oneoffprojects.walkers.IndelCountCovariates;

import net.sf.samtools.SAMRecord;
//g import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
//g import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: Jan 17, 2011
 * Time: 2:53:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class IndelPositionCovariate implements Covariate {

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
    }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {
        int cycle = offset;
        if( read.getReadNegativeStrandFlag() ) {
            cycle = read.getReadLength() - (offset + 1);
        }
        return cycle;
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }

    public void getValues(SAMRecord read, Comparable[] comparable) {
        for(int iii = 0; iii < read.getReadLength(); iii++) {
            comparable[iii] = getValue(read, iii); // BUGBUG: this can be optimized
        }
    }
}