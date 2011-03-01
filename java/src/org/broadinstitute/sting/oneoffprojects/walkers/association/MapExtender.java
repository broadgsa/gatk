package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.Map;

/**
 * @Author chartl
 * @Date 2011-02-23
 * Holds multiple map contexts for use in the regional association walker
 */
public class MapExtender {

    private MapHolder previous = null;
    private MapHolder current = null;

    public MapExtender() {
        // no need to do anything
    }

    public void set(MapHolder holder) {
        previous = current;
        current = holder;
    }

    public Map<Sample,StratifiedAlignmentContext> getPreviousContext() {
        return previous != null ? previous.getContext() : null;
    }

    public ReferenceContext getPreviousRef() {
        return previous != null ? previous.getRef() : null;
    }

    public RefMetaDataTracker getPreviousTracker() {
        return previous != null ? previous.getTracker() : null;
    }

    public Map<Sample,StratifiedAlignmentContext> getContext() {
        return current != null ? current.getContext() : null;
    }

    public ReferenceContext getReferenceContext() {
        return current != null ? current.getRef() : null;
    }

    public RefMetaDataTracker getTracker() {
        return current != null ? current.getTracker() : null;
    }
}
