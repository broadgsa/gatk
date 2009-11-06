package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 30, 2009
 */
public class ReadGroupCovariate implements Covariate{

    public ReadGroupCovariate() { // empty constructor is required by CovariateCounterWalker
    }

    public Comparable<?> getValue(SAMRecord read, int offset, char[] refBases) {
        return read.getReadGroup().getReadGroupId();
    }
    
    public Comparable<?> getValue(String str) {
    	return str;
    }

    public String toString() {
        return "Read Group";
    }
}


