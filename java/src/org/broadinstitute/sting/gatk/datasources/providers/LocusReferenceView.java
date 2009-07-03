package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.SAMSequenceRecord;
/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
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
     * Allow the user to pull reference info from any arbitrary region of the reference.
     * Assume the user has already performed all necessary bounds checking.
     * TODO: This function is nearly identical to that in the ReadReferenceView.  Merge the common functionality.
     * @param genomeLoc The locus.
     * @return A list of the bases starting at the start of the locus (inclusive) and ending
     *         at the end of the locus (inclusive).
     */
    public char[] getReferenceBases( GenomeLoc genomeLoc ) {
        SAMSequenceRecord sequenceInfo = reference.getSequenceDictionary().getSequence(genomeLoc.getContig());
        long stop = Math.min( genomeLoc.getStop(), sequenceInfo.getSequenceLength() );
        ReferenceSequence subsequence = reference.getSubsequenceAt(genomeLoc.getContig(),genomeLoc.getStart(),stop);
        return (StringUtil.bytesToString(subsequence.getBases()) + Utils.dupString('X', (int)(genomeLoc.getStop() - stop)) ).toCharArray();
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
