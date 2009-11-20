package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.alignment.Aligner;
import org.broadinstitute.sting.alignment.Alignment;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;

import java.util.Iterator;

/**
 * Align reads using BWA.
 *
 * @author mhanna
 * @version 0.1
 */
public abstract class BWAAligner implements Aligner {

    /**
     * Create a new BWAAligner.  Purpose of this call is to ensure that all BWA constructors accept the correct
     * parameters.
     * @param bwtFiles The many files representing BWTs persisted to disk.
     * @param configuration Configuration parameters for the alignment.
     */
    public BWAAligner(BWTFiles bwtFiles, BWAConfiguration configuration) { }

    /**
     * Update the configuration passed to the BWA aligner.
     * @param configuration New configuration to set.
     */    
    public abstract void updateConfiguration(BWAConfiguration configuration);
}
