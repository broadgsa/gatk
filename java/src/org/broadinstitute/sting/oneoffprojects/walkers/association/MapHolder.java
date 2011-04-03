package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.Map;

public class MapHolder {
    private RefMetaDataTracker tracker;
    private ReferenceContext ref;
    private Map<Sample, AlignmentContext> alignments;

    public MapHolder(RefMetaDataTracker t, ReferenceContext r, AlignmentContext a) {
        tracker = t;
        ref = r;
        alignments = StratifiedAlignmentContext.splitContextBySample(a);
    }

    public Map<Sample, AlignmentContext> getContext() {
        return alignments;
    }

    public ReferenceContext getRef() {
        return ref;
    }

    public RefMetaDataTracker getTracker() {
        return tracker;
    }
}