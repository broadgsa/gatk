package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.QualityUtils;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 */
public class QualityScoreCovariate implements Covariate {

    public final static String ORIGINAL_QUAL_ATTRIBUTE_TAG = "OQ";
    protected boolean USE_ORIGINAL_QUALS;

    public QualityScoreCovariate() { // empty constructor is required by CovariateCounterWalker
        USE_ORIGINAL_QUALS = false;
    }

    public QualityScoreCovariate(boolean originalQuals) {
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

        return quals[offset];
    }
    
    public Comparable<?> getValue(String str) {
        return Integer.parseInt( str );
    }
    
    public String toString() {
        return "Reported Quality Score";
    }
}
