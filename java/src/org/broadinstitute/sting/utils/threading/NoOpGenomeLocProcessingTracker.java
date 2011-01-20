package org.broadinstitute.sting.utils.threading;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Base class, and null tracker.  Always says that a GenomeLoc is ready for processing
 */
public class NoOpGenomeLocProcessingTracker extends GenomeLocProcessingTracker {
    protected NoOpGenomeLocProcessingTracker() {
        super(new ClosableReentrantLock());          // todo -- should be lighter weight
    }

    @Override
    public ProcessingLoc claimOwnership(GenomeLoc loc, String myName) {
        return new ProcessingLoc(loc, myName);
    }

    @Override
    protected List<ProcessingLoc> getProcessingLocs() {
        return Collections.emptyList();
    }

    @Override
    protected void registerNewLoc(ProcessingLoc loc) {
        ;
    }

    @Override
    protected List<ProcessingLoc> readNewLocs() {
        return Collections.emptyList();
    }
}
