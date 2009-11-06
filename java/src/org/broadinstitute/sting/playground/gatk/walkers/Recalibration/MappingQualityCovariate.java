package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 4, 2009
 */
public class MappingQualityCovariate implements Covariate {

    public MappingQualityCovariate() { // empty constructor is required by CovariateCounterWalker
    }

    public Comparable<?> getValue(SAMRecord read, int offset, char[] refBases) {
        return read.getMappingQuality();
    }
    
    public Comparable<?> getValue(String str) {
        return Integer.parseInt( str );
    }
    
    public String toString() {
        return "Mapping Quality Score";
    }
}
