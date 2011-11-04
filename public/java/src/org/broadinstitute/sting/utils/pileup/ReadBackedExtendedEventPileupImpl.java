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
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

public class ReadBackedExtendedEventPileupImpl extends AbstractReadBackedPileup<ReadBackedExtendedEventPileupImpl,ExtendedEventPileupElement> implements ReadBackedExtendedEventPileup {
    private int nInsertions;
    private int maxDeletionLength;      // cached value of the length of the longest deletion observed at the site

    public ReadBackedExtendedEventPileupImpl(GenomeLoc loc, List<ExtendedEventPileupElement> pileupElements) {
        super(loc,pileupElements);
    }

    public ReadBackedExtendedEventPileupImpl(GenomeLoc loc, PileupElementTracker<ExtendedEventPileupElement> tracker) {
        super(loc,tracker);
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc
     * @param pileup
     */
    public ReadBackedExtendedEventPileupImpl(GenomeLoc loc, List<ExtendedEventPileupElement> pileup, int size,
                                         int maxDeletionLength, int nInsertions, int nDeletions, int nMQ0Reads) {
        super(loc,pileup,size,nDeletions,nMQ0Reads);
        this.maxDeletionLength = maxDeletionLength;
        this.nInsertions = nInsertions;
    }

    // this is the good new one
    public ReadBackedExtendedEventPileupImpl(GenomeLoc loc, Map<String,? extends ReadBackedExtendedEventPileupImpl> pileupElementsBySample) {
        super(loc,pileupElementsBySample);
    }

    /**
     * Calculate cached sizes, nDeletion, and base counts for the pileup.  This calculation is done upfront,
     * so you pay the cost at the start, but it's more efficient to do this rather than pay the cost of calling
     * sizes, nDeletion, etc. over and over potentially.
     */
    @Override
    protected void calculateCachedData() {
        super.calculateCachedData();

        nInsertions = 0;
        nMQ0Reads = 0;

        for ( ExtendedEventPileupElement p : this.toExtendedIterable() ) {

            if ( p.isDeletion() ) {
                maxDeletionLength = Math.max(maxDeletionLength, p.getEventLength());
            } else {
                if ( p.isInsertion() ) nInsertions++;
            }
        }
    }

    @Override
    protected void addPileupToCumulativeStats(AbstractReadBackedPileup<ReadBackedExtendedEventPileupImpl,ExtendedEventPileupElement> pileup) {
        super.addPileupToCumulativeStats(pileup);
        ReadBackedExtendedEventPileup extendedEventPileup = ((ReadBackedExtendedEventPileup)pileup);
        this.nInsertions += extendedEventPileup.getNumberOfInsertions();
        this.maxDeletionLength += extendedEventPileup.getMaxDeletionLength();
    }

    @Override
    protected ReadBackedExtendedEventPileupImpl createNewPileup(GenomeLoc loc, PileupElementTracker<ExtendedEventPileupElement> tracker) {
        return new ReadBackedExtendedEventPileupImpl(loc,tracker);
    }

    @Override
    protected ExtendedEventPileupElement createNewPileupElement(GATKSAMRecord read, int offset) {
        throw new UnsupportedOperationException("Not enough information provided to create a new pileup element");
    }


    /**
     * Returns the number of insertion events in this pileup
     *
     * @return
     */
    @Override
    public int getNumberOfInsertions() {
        return nInsertions;
    }

    /** Returns the length of the longest deletion observed at the site this
     * pileup is associated with (NOTE: by convention, both insertions and deletions
     * are associated with genomic location immediately before the actual event). If
     * there are no deletions at the site, returns 0.
     * @return
     */
    @Override
    public int getMaxDeletionLength() {
        return maxDeletionLength;
    }

    public Iterable<ExtendedEventPileupElement> toExtendedIterable() {
        return new Iterable<ExtendedEventPileupElement>() {
            public Iterator<ExtendedEventPileupElement> iterator() { return pileupElementTracker.iterator(); }
        };
    }

    /**
     * Returns an array of the events in this pileup ('I', 'D', or '.'). Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    @Override
    public byte[] getEvents() {
        byte[] v = new byte[getNumberOfElements()];
        int i = 0;
        for ( ExtendedEventPileupElement e : this.toExtendedIterable() ) {
            switch ( e.getType() ) {
                case INSERTION: v[i] = 'I'; break;
                case DELETION: v[i] = 'D'; break;
                case NOEVENT: v[i] = '.'; break;
                default: throw new ReviewedStingException("Unknown event type encountered: "+e.getType());
            }
            i++;
        }
        return v;
    }    

    /** A shortcut for getEventStringsWithCounts(null);
     *
     * @return
     */
    @Override
    public List<Pair<String,Integer>> getEventStringsWithCounts() {
        return getEventStringsWithCounts(null);
    }

    @Override
    public String getShortPileupString() {
        // In the pileup format, each extended event line has genomic position (chromosome name and offset),
        // reference "base" (always set to "E" for E(xtended)), then 'I','D', or '.' symbol for each read representing
        // insertion, deletion or no-event, respectively.
        return String.format("%s %s E %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                new String(getEvents()) );
    }

    /** Returns String representation of all distinct extended events (indels) at the site along with
     * observation counts (numbers of reads) for each distinct event. If refBases is null, a simple string representation for
     * deletions will be generated as "<length>D" (i.e. "5D"); if the reference bases are provided, the actual
     * deleted sequence will be used in the string representation (e.g. "-AAC").
      * @param refBases reference bases, starting with the current locus (i.e. the one immediately before the indel), and
     * extending far enough to accomodate the longest deletion (i.e. size of refBases must be at least 1+<length of longest deletion>)
     * @return list of distinct events; first element of a pair is a string representation of the event, second element
     * gives the number of reads, in which that event was observed
     */
    @Override
    public List<Pair<String,Integer>> getEventStringsWithCounts(byte[] refBases) {
        Map<String, Integer> events = new HashMap<String,Integer>();

        for ( ExtendedEventPileupElement e : this.toExtendedIterable() ) {
            Integer cnt;
            String indel = null;
            switch ( e.getType() ) {
                case INSERTION:
                    indel = "+"+e.getEventBases();
                    break;
                case DELETION:
                    indel = getDeletionString(e.getEventLength(),refBases);
                    break;
                case NOEVENT: continue;
                default: throw new ReviewedStingException("Unknown event type encountered: "+e.getType());
            }

            cnt = events.get(indel);
            if ( cnt == null ) events.put(indel,1);
            else events.put(indel,cnt.intValue()+1);
        }

        List<Pair<String,Integer>> eventList = new ArrayList<Pair<String,Integer>>(events.size());
        for ( Map.Entry<String,Integer> m : events.entrySet() ) {
            eventList.add( new Pair<String,Integer>(m.getKey(),m.getValue()));
        }
        return eventList;
    }

    /**
     * Builds string representation of the deletion event. If refBases is null, the representation will be
     * "<length>D" (e.g. "5D"); if the reference bases are provided, a verbose representation (e.g. "-AAC")
     *  will be generated. NOTE: refBases must start with the base prior to the actual deletion (i.e. deleted
     * base(s) are refBase[1], refBase[2], ...), and the length of the passed array must be sufficient to accomodate the
     * deletion length (i.e. size of refBase must be at least length+1).
     * @param length
     * @param refBases
     * @return
     */
    private String getDeletionString(int length, byte[] refBases) {
        if ( refBases == null ) {
            return  Integer.toString(length)+"D"; // if we do not have reference bases, we can only report something like "5D"
        } else {
            return "-"+new String(refBases,1,length).toUpperCase();
        }
    }
}