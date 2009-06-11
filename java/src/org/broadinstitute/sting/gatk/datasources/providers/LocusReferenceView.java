package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;
/**
 * User: hanna
 * Date: May 22, 2009
 * Time: 12:24:23 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Provides access to the portion of the reference covering a single locus.
 */
public class LocusReferenceView extends ReferenceView {
    /**
     * Bound the reference view to make sure all accesses are within the shard.
     */
    private final GenomeLoc bounds;

    /**
     * Track the reference sequence and the last point accessed.  Used to
     * track state when traversing over the reference.
     */
    private ReferenceSequence referenceSequence;

    /**
     * Create a new locus reference view.
     * @param provider source for locus data.
     */
    public LocusReferenceView( ShardDataProvider provider ) {
        super( provider );
        bounds = provider.getShard().getGenomeLoc();
        this.referenceSequence = reference.getSubsequenceAt( bounds.getContig(),
                                                             bounds.getStart(),
                                                             bounds.getStop() );        
    }

    /**
     * Gets the reference base associated with this particular point on the genome.
     * @param genomeLoc Region for which to retrieve the base.  GenomeLoc must represent a 1-base region.
     * @return The base at the position represented by this genomeLoc.
     */
    public char getReferenceBase( GenomeLoc genomeLoc ) {
        validateLocation( genomeLoc );
        int offset = (int)(genomeLoc.getStart() - bounds.getStart());
        return StringUtil.bytesToString( referenceSequence.getBases(), offset, 1 ).charAt(0);
    }

    /**
     * Validates that the genomeLoc is one base wide and is in the reference sequence.
     * @param genomeLoc location to verify.
     */
    private void validateLocation( GenomeLoc genomeLoc ) throws InvalidPositionException {
        //
        if( !genomeLoc.isSingleBP() )
            throw new InvalidPositionException(
                    String.format("Requested position larger than one base; start = %d, stop = %d", genomeLoc.getStart(), genomeLoc.getStop()));
        if( !bounds.containsP(genomeLoc) )
            throw new InvalidPositionException(
                    String.format("Requested position %s not within interval %s", genomeLoc, bounds));
    }
}
