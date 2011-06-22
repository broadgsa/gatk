package org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/8/11
 * Time: 11:59 AM
 * To change this template use File | Settings | File Templates.
 */
public class ProperPair extends BinaryFeatureAggregator {

    public ProperPair(RFAArgumentCollection collection) {
        super(collection);
    }

    public Boolean extractFeature(SAMRecord record) {
        return record.getProperPairFlag();
    }

    public boolean featureDefined(SAMRecord record) {
        return record.getReadPairedFlag();
    }
}
