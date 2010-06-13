/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMRecord;

/**
 * A clean interface for working with extended event pileups.
 *
 * @author mhanna
 * @version 0.1
 */
public interface ReadBackedExtendedEventPileup extends Iterable<ExtendedEventPileupElement> {
    /**
     * Returns a new ReadBackedPileup that is free of mapping quality zero reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no MQ0 reads in the pileup.
     *
     * @return
     */
    public ReadBackedExtendedEventPileup getPileupWithoutMappingQualityZeroReads();
    
    public ReadBackedExtendedEventPileup getMappingFilteredPileup( int minMapQ );

    /**
     * Returns a pileup randomly downsampled to the desiredCoverage.
     *
     * @param desiredCoverage
     * @return
     */
    public ReadBackedExtendedEventPileup getDownsampledPileup(int desiredCoverage);

    /**
     * Returns the number of deletion events in this pileup
     *
     * @return
     */
    public int getNumberOfDeletions();

    /**
     * Returns the number of insertion events in this pileup
     *
     * @return
     */
    public int getNumberOfInsertions();

    /** Returns the length of the longest deletion observed at the site this
     * pileup is associated with (NOTE: by convention, both insertions and deletions
     * are associated with genomic location immediately before the actual event). If
     * there are no deletions at the site, returns 0.
     * @return
     */
    public int getMaxDeletionLength();

    /**
     * Returns the number of mapping quality zero reads in this pileup.
     * @return
     */
    public int getNumberOfMappingQualityZeroReads();

    /**
     * @return the number of elements in this pileup
     */
    public int size();

    /**
     * @return the location of this pileup
     */
    public GenomeLoc getLocation();

    public String getShortPileupString();

    /**
     * Returns a list of the reads in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    public List<SAMRecord> getReads();

    /**
     * Returns a list of the offsets in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    public List<Integer> getOffsets();

    /**
     * Returns an array of the events in this pileup ('I', 'D', or '.'). Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    public byte[] getEvents();

    /** A shortcut for getEventStringsWithCounts(null);
     *
     * @return
     */
    public List<Pair<String,Integer>> getEventStringsWithCounts();

    /** Returns String representation of all distinct extended events (indels) at the site along with
     * observation counts (numbers of reads) for each distinct event. If refBases is null, a simple string representation for
     * deletions will be generated as "<length>D" (i.e. "5D"); if the reference bases are provided, the actual
     * deleted sequence will be used in the string representation (e.g. "-AAC").
      * @param refBases reference bases, starting with the current locus (i.e. the one immediately before the indel), and
     * extending far enough to accomodate the longest deletion (i.e. size of refBases must be at least 1+<length of longest deletion>)
     * @return list of distinct events; first element of a pair is a string representation of the event, second element
     * gives the number of reads, in which that event was observed
     */
    public List<Pair<String,Integer>> getEventStringsWithCounts(byte[] refBases);

    /**
     * Get an array of the mapping qualities
     * @return
     */
    public byte[] getMappingQuals();    
}
