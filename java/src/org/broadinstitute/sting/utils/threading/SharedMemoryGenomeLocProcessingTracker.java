package org.broadinstitute.sting.utils.threading;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.List;

/**
 * For algorithmic testing purposes only.  Uses synchronization to keep a consistent
 * processing list in shared memory.
 */
public class SharedMemoryGenomeLocProcessingTracker extends GenomeLocProcessingTracker {
    List<ProcessingLoc> processingLocs = new ArrayList<ProcessingLoc>();

    public synchronized ProcessingLoc claimOwnership(GenomeLoc loc, String myName) {
        // processingLocs is a shared memory synchronized object, and this
        // method is synchonized, so we can just do our processing
        ProcessingLoc owner = super.findOwner(loc);

        if ( owner == null ) { // we are unowned
            owner = new ProcessingLoc(loc, myName);
            processingLocs.add(owner);
        }

        return owner;
    }

    protected synchronized List<ProcessingLoc> getProcessingLocs() {
        return processingLocs;
    }
}
