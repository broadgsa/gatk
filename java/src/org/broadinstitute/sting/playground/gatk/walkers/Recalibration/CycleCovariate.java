package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 30, 2009
 */
public class CycleCovariate implements Covariate {

    public String platform;

    public CycleCovariate() { // empty constructor is required by CovariateCounterWalker
        platform = null;
    }

    public CycleCovariate(String _platform) {
        platform = _platform;
    }

    public Comparable<?> getValue(SAMRecord read, int offset, char[] refBases) {
        //BUGBUG: assumes Solexia platform
        int cycle = offset;
        if( read.getReadNegativeStrandFlag() ) {
            cycle = read.getReadLength() - (offset + 1);
        }
        return cycle;
    }
    
    public Comparable<?> getValue(String str) {
        return Integer.parseInt( str );
    }

    public String toString() {
        return "Cycle";
    }
}