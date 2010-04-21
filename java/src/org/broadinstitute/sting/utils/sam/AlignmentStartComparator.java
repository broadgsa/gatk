package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMRecord;

import java.util.Comparator;

/**
 * Compares two SAMRecords only the basis on alignment start.  Note that
 * comparisons are performed ONLY on the basis of alignment start; any
 * two SAM records with the same alignment start will be considered equal.
 *
 * Unmapped alignments will all be considered equal.
 *
 * @author mhanna
 * @version 0.1
 */
public class AlignmentStartComparator implements Comparator<SAMRecord> {
    public int compare(SAMRecord lhs, SAMRecord rhs) {
        if(!lhs.getReferenceIndex().equals(rhs.getReferenceIndex()))
            return lhs.getReferenceIndex() - rhs.getReferenceIndex();

        // Note: no integer overflow here because alignment starts are >= 0.
        return lhs.getAlignmentStart() - rhs.getAlignmentStart();
    }
}
