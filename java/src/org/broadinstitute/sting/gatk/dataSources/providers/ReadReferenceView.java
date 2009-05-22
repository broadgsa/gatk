package org.broadinstitute.sting.gatk.dataSources.providers;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.StringUtil;
import edu.mit.broad.picard.reference.ReferenceSequence;
/**
 * User: hanna
 * Date: May 22, 2009
 * Time: 12:36:14 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Provides access to the reference over a single read.
 */

public class ReadReferenceView extends ReferenceView {
    /**
     * Create a view of the reference with respect to a single read.
     * @param provider 
     */
    public ReadReferenceView( ShardDataProvider provider ) {
        super( provider );
    }

    /**
     * Gets the bases of the reference that are aligned to the given read.
     * @param read the read for which to extract reference information.
     * @return The bases corresponding to this read, or null if the read is unmapped.
     *         If the alignment goes off the end of the contig, return just the portion
     *         mapped to the reference.
     */
    public char[] getReferenceBases( SAMRecord read ) {
        if( read.getReadUnmappedFlag() )
            return null;

        String contig = read.getReferenceName();
        int start = read.getAlignmentStart();
        int stop = read.getAlignmentEnd();

        SAMSequenceRecord sequenceRecord = reference.getSequenceDictionary().getSequence(contig);
        if( stop > sequenceRecord.getSequenceLength() )
            stop = sequenceRecord.getSequenceLength();

        ReferenceSequence alignmentToReference = reference.getSubsequenceAt( contig, start, stop );
        return StringUtil.bytesToString(alignmentToReference.getBases()).toCharArray();
    }

}
