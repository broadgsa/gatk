package org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 1:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class AberrantInsertSize extends BinaryFeatureAggregator {

    private int min;
    private int max;

    public AberrantInsertSize(RFAArgumentCollection col) {
        super(col);
        min = col.lowInsertSize;
        max = col.highInsertSize;
    }

    public Boolean extractFeature(SAMRecord rec) {
        return Math.abs(rec.getInferredInsertSize()) > max || Math.abs(rec.getInferredInsertSize()) < min;
    }

    public boolean featureDefined(SAMRecord rec) {
        return rec.getReadPairedFlag() && rec.getProperPairFlag();
    }
}
