package org.broadinstitute.sting.utils.threading;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Base class, and null tracker.  Always says that a GenomeLoc is ready for processing.  It is
 * critical that this class already return that a loc is owned, no matter if it's been seen before,
 * etc.  ReadShards can differ in their contents but have the same "unmapped" genome loc
 */
public class NoOpGenomeLocProcessingTracker extends GenomeLocProcessingTracker {
    public NoOpGenomeLocProcessingTracker() {
        super(new ClosableReentrantLock(), null);
    }

    @Override
    protected void registerNewLocs(Collection<ProcessingLoc> loc) {
        ;
    }

    @Override
    protected List<ProcessingLoc> readNewLocs() {
        return Collections.emptyList();
    }
}
