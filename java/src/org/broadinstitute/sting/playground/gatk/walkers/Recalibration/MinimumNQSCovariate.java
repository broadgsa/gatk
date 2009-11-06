package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.QualityUtils;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 4, 2009
 */
public class MinimumNQSCovariate implements Covariate {

    public final static String ORIGINAL_QUAL_ATTRIBUTE_TAG = "OQ";
    protected boolean USE_ORIGINAL_QUALS;

    public MinimumNQSCovariate() { // empty constructor is required by CovariateCounterWalker
        USE_ORIGINAL_QUALS = false;
    }

    public MinimumNQSCovariate(boolean originalQuals) {
        USE_ORIGINAL_QUALS = originalQuals;
    }

    public Comparable<?> getValue(SAMRecord read, int offset, char[] refBases) {
        byte[] quals = read.getBaseQualities();
        if ( USE_ORIGINAL_QUALS && read.getAttribute(ORIGINAL_QUAL_ATTRIBUTE_TAG) != null ) {
            Object obj = read.getAttribute(ORIGINAL_QUAL_ATTRIBUTE_TAG);
            if ( obj instanceof String )
                quals = QualityUtils.fastqToPhred((String)obj);
            else {
                throw new RuntimeException(String.format("Value encoded by %s in %s isn't a string!", ORIGINAL_QUAL_ATTRIBUTE_TAG, read.getReadName()));
            }
        }

        int minQual = quals[0];
        for ( int qual : quals ) {
            if( qual < minQual ) {
                minQual = qual;
            }
        }
        return minQual;
    }
    
    public Comparable<?> getValue(String str) {
        return Integer.parseInt( str );
    }

    public String toString() {
        return "Minimum Neighborhood Quality Score";
    }
}
