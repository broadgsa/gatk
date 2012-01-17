package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 12/7/11
 */

@By(DataSource.READS)
@Requires({DataSource.READS, DataSource.REFERENCE_BASES})
@PartitionBy(PartitionType.READ)
public abstract class ActiveRegionWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    // Determine active status over the AlignmentContext
    public abstract boolean isActive(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context);

    // Map over the ActiveRegion
    public abstract MapType map(final ActiveRegion activeRegion, final ReadMetaDataTracker metaDataTracker);
}
