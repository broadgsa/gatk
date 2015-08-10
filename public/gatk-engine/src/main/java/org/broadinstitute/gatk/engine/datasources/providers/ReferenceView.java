/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
/**
 * User: hanna
 * Date: May 22, 2009
 * Time: 12:19:17 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A view into the reference backing this shard.
 */
public class ReferenceView implements View {
    /**
     * The parser, used to create and parse GenomeLocs.
     */
    protected final GenomeLocParser genomeLocParser;

    /**
     * The source of reference data.
     */
    protected IndexedFastaSequenceFile reference = null;

    /**
     * Create a new ReferenceView.
     * @param provider
     */
    public ReferenceView( ShardDataProvider provider ) {
        this.genomeLocParser = provider.getGenomeLocParser();
        this.reference = provider.getReference();
    }

    /**
     * Reference views don't conflict with anything else.
     * @return Empty list.
     */
    public Collection<Class<? extends View>> getConflictingViews() { return Collections.emptyList(); }

    /**
     * Deinitialize pointers for fast fail.  Someone else will handle file management.
     */
    public void close() {
        reference = null;
    }

    /**
     * Allow the user to pull reference info from any arbitrary region of the reference.
     * If parts of the reference don't exist, mark them in the char array with 'X'es.
     * @param genomeLoc The locus.
     * @return A list of the bases starting at the start of the locus (inclusive) and ending
     *         at the end of the locus (inclusive).
     */
    final static int BUFFER = 10000;
    final static byte[] Xs = new byte[BUFFER];
    static {
        Arrays.fill(Xs, (byte)'X');
    }

    protected byte[] getReferenceBases( SAMRecord read ) {
        return getReferenceBases(genomeLocParser.createGenomeLoc(read));

    }

    protected byte[] getReferenceBases( GenomeLoc genomeLoc ) {
        SAMSequenceRecord sequenceInfo = reference.getSequenceDictionary().getSequence(genomeLoc.getContig());

        long start = genomeLoc.getStart();
        long stop = Math.min( genomeLoc.getStop(), sequenceInfo.getSequenceLength() );

        // Read with no aligned bases?  Return an empty array.
        if(stop - start + 1 == 0)
            return new byte[0];

        ReferenceSequence subsequence = reference.getSubsequenceAt(genomeLoc.getContig(), start, stop);

        int overhang = (int)(genomeLoc.getStop() - stop);
        if ( overhang > 0 ) {
            if ( overhang > BUFFER ) // todo -- this is a bit dangerous
                throw new ReviewedGATKException("Insufficient buffer size for Xs overhanging genome -- expand BUFFER");
            byte[] all = new byte[subsequence.getBases().length + overhang];
            System.arraycopy(subsequence.getBases(), 0, all, 0, subsequence.getBases().length);
            System.arraycopy(Xs, 0, all, subsequence.getBases().length, overhang);
            return all;
        } else {
            // fast path
            return subsequence.getBases();
        }
    }
}
