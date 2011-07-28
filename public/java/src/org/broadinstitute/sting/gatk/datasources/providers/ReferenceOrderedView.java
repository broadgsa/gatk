package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;

public interface ReferenceOrderedView extends View {
    RefMetaDataTracker getReferenceOrderedDataAtLocus( GenomeLoc loc, ReferenceContext refContext );
}
