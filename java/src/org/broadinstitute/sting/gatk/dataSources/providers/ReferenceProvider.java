package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 8, 2009
 * Time: 5:01:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceProvider {
    private IndexedFastaSequenceFile sequenceFile;
    private Shard shard;

    /**
     * Track the reference sequence and the last point accessed.  Used to
     * track state when traversing over the reference.
     */
    private ReferenceSequence referenceSequence;
    private GenomeLoc referenceInterval;

    /**
     * Create a new reference provider supplying data from the given reference.
     * @param sequenceFile Reference file to use.
     * @param shard Shard over which to retrieve data.
     */
    public ReferenceProvider( IndexedFastaSequenceFile sequenceFile, Shard shard ) {
        this.sequenceFile = sequenceFile;
        this.shard = shard;
    }

    /**
     * Gets the reference base at a single point.
     * @param genomeLoc The location at which to fetch the reference base.
     * @return The character representing the reference base.
     * @throws InvalidPositionException in case the position is invalid.
     */
    public char getReferenceBase( GenomeLoc genomeLoc ) throws InvalidPositionException {
        if( referenceSequence == null )
            lazyInitializeLocusAccess();            

        validateLocation( genomeLoc );
        int offset = (int)(genomeLoc.getStart() - referenceInterval.getStart());
        return StringUtil.bytesToString( referenceSequence.getBases(), offset, 1 ).charAt(0);
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

        SAMSequenceRecord sequenceRecord = sequenceFile.getSequenceDictionary().getSequence(contig);
        if( stop > sequenceRecord.getSequenceLength() )
            stop = sequenceRecord.getSequenceLength();

        ReferenceSequence alignmentToReference = sequenceFile.getSubsequenceAt( contig, start, stop );
        return StringUtil.bytesToString(alignmentToReference.getBases()).toCharArray();
    }

    /**
     * Perform a lazy initialization of access to the locus.  Sets up the reference sequence and
     * limits the user to work only at that sequence.
     */
    private void lazyInitializeLocusAccess() {
        GenomeLoc position = shard.getGenomeLoc();
        this.referenceSequence = sequenceFile.getSubsequenceAt( position.getContig(),
                                                                position.getStart(),
                                                                position.getStop() );
        this.referenceInterval = position;
    }

    /**
     * Validates that the genomeLoc is one base wide and is in the reference sequence.
     * @param genomeLoc location to verify.
     */
    private void validateLocation( GenomeLoc genomeLoc ) throws InvalidPositionException {
        //
        if( !referenceInterval.containsP(genomeLoc) )
            throw new InvalidPositionException(
                    String.format("Requested position %s not within interval %s", genomeLoc, referenceInterval));
        if( genomeLoc.getStart() != genomeLoc.getStop() )
            throw new InvalidPositionException(
                    String.format("Requested position larger than one base; start = %d, stop = %d", genomeLoc.getStart(), genomeLoc.getStop()));
    }
}
