package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.HashMap;
import java.util.Map;

public class MapHolder {
    private RefMetaDataTracker tracker;
    private ReferenceContext ref;
    private Map<String, AlignmentContext> alignments;

    public MapHolder(RefMetaDataTracker t, ReferenceContext r, AlignmentContext a) {
        tracker = t;
        ref = r;
        alignments = AlignmentContextUtils.splitContextBySampleName(a);
    }

    public Map<Sample, AlignmentContext> getContext(Map<String,Sample> stringSampleMap) {
        Map<Sample,AlignmentContext> mappedContexts = new HashMap<Sample,AlignmentContext>(alignments.size());
        for ( Map.Entry<String,AlignmentContext> entry : alignments.entrySet() ) {
            mappedContexts.put(stringSampleMap.get(entry.getKey()),entry.getValue());
        }
        return mappedContexts;
    }

    public ReferenceContext getRef() {
        return ref;
    }

    public RefMetaDataTracker getTracker() {
        return tracker;
    }
}