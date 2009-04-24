package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 8, 2009
 * Time: 5:01:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceProvider {
    private ReferenceSequence referenceSequence;
    private GenomeLoc referenceInterval;

    public ReferenceProvider( IndexedFastaSequenceFile sequenceFile, GenomeLoc position ) {
        this.referenceSequence = sequenceFile.getSubsequenceAt( position.getContig(),
                                                                position.getStart(),
                                                                position.getStop() );
        this.referenceInterval = position;
    }

    public char getReferenceBase( GenomeLoc genomeLoc ) throws InvalidPositionException {
        validateLocation( genomeLoc );
        int offset = (int)(genomeLoc.getStart() - referenceInterval.getStart());
        return StringUtil.bytesToString( referenceSequence.getBases(), offset, 1 ).charAt(0);
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
