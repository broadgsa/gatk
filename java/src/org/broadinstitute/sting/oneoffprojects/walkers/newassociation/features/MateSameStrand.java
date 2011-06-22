package org.broadinstitute.sting.oneoffprojects.walkers.newassociation.features;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.RFAArgumentCollection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/4/11
 * Time: 1:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateSameStrand extends BinaryFeatureAggregator {

    public Boolean extractFeature(SAMRecord record) {
        return record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
    }

    public boolean featureDefined(SAMRecord record) {
        return record.getReadPairedFlag() && ! record.getMateUnmappedFlag() && record.getReferenceIndex().equals(record.getMateReferenceIndex());
    }


    public MateSameStrand(RFAArgumentCollection col) {
        super(col);
    }
}
