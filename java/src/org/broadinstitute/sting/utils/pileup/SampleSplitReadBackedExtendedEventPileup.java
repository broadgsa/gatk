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
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.*;

import net.sf.samtools.SAMRecord;

/**
 * Extended event pileups, split out by sample.
 *
 * @author mhanna
 * @version 0.1
 */
public class SampleSplitReadBackedExtendedEventPileup implements ReadBackedExtendedEventPileup {
    private GenomeLoc loc = null;
    private Map<String,ReadBackedExtendedEventPileup> pileupBySample = null;

    private int size = 0;                   // cached value of the size of the pileup
    private int maxDeletionLength = 0;      // cached value of the length of the longest deletion observed at the site
    private int nDeletions = 0;             // cached value of the number of deletions
    private int nInsertions = 0;
    private int nMQ0Reads = 0;              // cached value of the number of MQ0 reads

    public SampleSplitReadBackedExtendedEventPileup(GenomeLoc loc) {
        this(loc,new HashMap<String,ReadBackedExtendedEventPileup>());
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc Location of this pileup.
     * @param pileup Pileup data, split by sample name.
     */
    public SampleSplitReadBackedExtendedEventPileup(GenomeLoc loc, Map<String,ReadBackedExtendedEventPileup> pileup) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in UnifiedReadBackedPileup");
        if ( pileup == null ) throw new StingException("Illegal null pileup in UnifiedReadBackedPileup");

        this.loc = loc;
        this.pileupBySample = pileup;
        for(ReadBackedExtendedEventPileup pileupBySample: pileup.values()) {
            this.size += pileupBySample.size();
            this.nDeletions += pileupBySample.getNumberOfDeletions();
            this.nMQ0Reads += pileupBySample.getNumberOfMappingQualityZeroReads();
        }
    }



    private void addPileup(final String sampleName, ReadBackedExtendedEventPileup pileup) {
        pileupBySample.put(sampleName,pileup);

        size += pileup.size();
        nDeletions += pileup.getNumberOfDeletions();
        nMQ0Reads += pileup.getNumberOfMappingQualityZeroReads();
    }

    @Override
    public ReadBackedExtendedEventPileup getPileupWithoutMappingQualityZeroReads() {
        if ( getNumberOfMappingQualityZeroReads() > 0 ) {
            SampleSplitReadBackedExtendedEventPileup filteredPileup = new SampleSplitReadBackedExtendedEventPileup(loc);

            for (Map.Entry<String,ReadBackedExtendedEventPileup> entry: pileupBySample.entrySet()) {
                String sampleName = entry.getKey();
                ReadBackedExtendedEventPileup pileup = entry.getValue();
                filteredPileup.addPileup(sampleName,pileup.getPileupWithoutMappingQualityZeroReads());
            }

            return filteredPileup;
        } else {
            return this;
        }
    }

    @Override
    public ReadBackedExtendedEventPileup getMappingFilteredPileup( int minMapQ ) {
        SampleSplitReadBackedExtendedEventPileup filteredPileup = new SampleSplitReadBackedExtendedEventPileup(loc);

        for (Map.Entry<String,ReadBackedExtendedEventPileup> entry: pileupBySample.entrySet()) {
            String sampleName = entry.getKey();
            ReadBackedExtendedEventPileup pileup = entry.getValue();
            filteredPileup.addPileup(sampleName,pileup.getMappingFilteredPileup(minMapQ));
        }

        return filteredPileup;
    }

    /**
     * Returns a pileup randomly downsampled to the desiredCoverage.
     *
     * @param desiredCoverage
     * @return
     */
    @Override
    public ReadBackedExtendedEventPileup getDownsampledPileup(int desiredCoverage) {
        if ( size() <= desiredCoverage )
            return this;

        // randomly choose numbers corresponding to positions in the reads list
        Random generator = new Random();
        TreeSet<Integer> positions = new TreeSet<Integer>();
        for ( int i = 0; i < desiredCoverage; /* no update */ ) {
            if ( positions.add(generator.nextInt(size)) )
                i++;
        }

        SampleSplitReadBackedExtendedEventPileup downsampledPileup = new SampleSplitReadBackedExtendedEventPileup(getLocation());
        int current = 0;

        for(Map.Entry<String,ReadBackedExtendedEventPileup> entry: this.pileupBySample.entrySet()) {
            String sampleName = entry.getKey();
            ReadBackedExtendedEventPileup pileup = entry.getValue();

            List<ExtendedEventPileupElement> filteredPileup = new ArrayList<ExtendedEventPileupElement>();
            for(ExtendedEventPileupElement p: pileup) {
                if(positions.contains(current))
                    filteredPileup.add(p);
            }

            if(!filteredPileup.isEmpty())
                downsampledPileup.addPileup(sampleName,new UnifiedReadBackedExtendedEventPileup(loc,filteredPileup));

            current++;
        }

        return downsampledPileup;
    }

    @Override
    public Iterator<ExtendedEventPileupElement> iterator() {
        return new ExtendedEventCastingIterator(new MergingPileupElementIterator(pileupBySample));
    }

    /**
     * Returns the number of deletion events in this pileup
     *
     * @return
     */
    @Override
    public int getNumberOfDeletions() {
        return nDeletions;
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
        int maxDeletionLength = 0;
        for(ReadBackedExtendedEventPileup pileup: pileupBySample.values())
            maxDeletionLength = Math.max(maxDeletionLength,pileup.getMaxDeletionLength());
        return maxDeletionLength;
    }

    /**
     * Returns the number of mapping quality zero reads in this pileup.
     * @return
     */
    @Override
    public int getNumberOfMappingQualityZeroReads() {
        return nMQ0Reads;
    }

    /**
     * @return the number of elements in this pileup
     */
    @Override
    public int size() {
        return size;
    }

    /**
     * @return the location of this pileup
     */
    @Override
    public GenomeLoc getLocation() {
        return loc;
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

    // --------------------------------------------------------
    //
    // Convenience functions that may be slow
    //
    // --------------------------------------------------------

    /**
     * Returns a list of the reads in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    @Override
    public List<SAMRecord> getReads() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>(size());
        for ( ExtendedEventPileupElement pile : this ) { reads.add(pile.getRead()); }
        return reads;
    }

    /**
     * Returns a list of the offsets in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    @Override
    public List<Integer> getOffsets() {
        List<Integer> offsets = new ArrayList<Integer>(size());
        for ( ExtendedEventPileupElement pile : this ) { offsets.add(pile.getOffset()); }
        return offsets;
    }

    /**
     * Returns an array of the events in this pileup ('I', 'D', or '.'). Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    @Override
    public byte[] getEvents() {
        byte[] v = new byte[size()];
        int i = 0;
        for ( ExtendedEventPileupElement e : this ) {
            switch ( e.getType() ) {
                case INSERTION: v[i] = 'I'; break;
                case DELETION: v[i] = 'D'; break;
                case NOEVENT: v[i] = '.'; break;
                default: throw new StingException("Unknown event type encountered: "+e.getType());
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

        for ( ExtendedEventPileupElement e : this ) {
            Integer cnt;
            String indel = null;
            switch ( e.getType() ) {
                case INSERTION:
                    indel = "+"+e.getEventBases();
                    break;
                case DELETION:
                    indel = UnifiedReadBackedExtendedEventPileup.getDeletionString(e.getEventLength(),refBases);
                    break;
                case NOEVENT: continue;
                default: throw new StingException("Unknown event type encountered: "+e.getType());
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
     * Get an array of the mapping qualities
     * @return
     */
    @Override
    public byte[] getMappingQuals() {
       byte[] v = new byte[size()];
       int i = 0;
       for ( ExtendedEventPileupElement e : this ) { v[i++] = (byte)e.getRead().getMappingQuality(); }
       return v;
    }    

    private class ExtendedEventCastingIterator implements Iterator<ExtendedEventPileupElement> {
        private final Iterator<PileupElement> wrappedIterator;

        public ExtendedEventCastingIterator(Iterator<PileupElement> iterator) {
            this.wrappedIterator = iterator;
        }

        public boolean hasNext() { return wrappedIterator.hasNext(); }
        public ExtendedEventPileupElement next() { return (ExtendedEventPileupElement)wrappedIterator.next(); }
        public void remove() { wrappedIterator.remove(); }
    }

}
