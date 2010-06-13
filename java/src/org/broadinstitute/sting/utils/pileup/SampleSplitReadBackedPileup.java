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

import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMRecord;

import java.util.*;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.iterators.IterableIterator;

/**
 * A read-backed pileup that internally divides the dataset based
 * on sample, for super efficient per-sample access.
 *
 * TODO: there are a few functions that are shared between UnifiedReadBackedPileup and SampleSplitReadBackedPileup.
 * TODO: refactor these into a common class. 
 *
 * @author mhanna
 * @version 0.1
 */
public class SampleSplitReadBackedPileup implements ReadBackedPileup {
    private GenomeLoc loc = null;
    private Map<String,ReadBackedPileup> pileupBySample = null;

    private int size = 0;                   // cached value of the size of the pileup
    private int nDeletions = 0;             // cached value of the number of deletions
    private int nMQ0Reads = 0;              // cached value of the number of MQ0 reads

    public SampleSplitReadBackedPileup(GenomeLoc loc) {
        this(loc,new HashMap<String,ReadBackedPileup>());
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc Location of this pileup.
     * @param pileup Pileup data, split by sample name.
     */
    public SampleSplitReadBackedPileup(GenomeLoc loc, Map<String,ReadBackedPileup> pileup) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in UnifiedReadBackedPileup");
        if ( pileup == null ) throw new StingException("Illegal null pileup in UnifiedReadBackedPileup");

        this.loc = loc;
        this.pileupBySample = pileup;
        for(ReadBackedPileup pileupBySample: pileup.values()) {
            this.size += pileupBySample.size();
            this.nDeletions += pileupBySample.getNumberOfDeletions();
            this.nMQ0Reads += pileupBySample.getNumberOfMappingQualityZeroReads();
        }
    }

    private void addPileup(final String sampleName, ReadBackedPileup pileup) {
        pileupBySample.put(sampleName,pileup);

        size += pileup.size();
        nDeletions += pileup.getNumberOfDeletions();
        nMQ0Reads += pileup.getNumberOfMappingQualityZeroReads();
    }

    /**
     * Returns a new ReadBackedPileup that is free of deletion spanning reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no deletions in the pileup.
     *
     * @return
     */
    @Override
    public ReadBackedPileup getPileupWithoutDeletions() {
        if ( getNumberOfDeletions() > 0 ) {
            SampleSplitReadBackedPileup filteredPileup = new SampleSplitReadBackedPileup(loc);

            for (Map.Entry<String,ReadBackedPileup> entry: pileupBySample.entrySet()) {
                String sampleName = entry.getKey();
                ReadBackedPileup pileup = entry.getValue();
                filteredPileup.addPileup(sampleName,pileup.getPileupWithoutDeletions());
            }

            return filteredPileup;
        } else {
            return this;
        }
    }

    /**
     * Returns a new ReadBackedPileup where only one read from an overlapping read
     * pair is retained.  If the two reads in question disagree to their basecall,
     * neither read is retained.  If they agree on the base, the read with the higher
     * quality observation is retained
     *
     * @return the newly filtered pileup
     */
    @Override
    public ReadBackedPileup getOverlappingFragmentFilteredPileup() {
        SampleSplitReadBackedPileup filteredPileup = new SampleSplitReadBackedPileup(loc);

        for (Map.Entry<String,ReadBackedPileup> entry: pileupBySample.entrySet()) {
            String sampleName = entry.getKey();
            ReadBackedPileup pileup = entry.getValue();
            filteredPileup.addPileup(sampleName,pileup.getOverlappingFragmentFilteredPileup());
        }

        return filteredPileup;
    }    

    /**
     * Returns a new ReadBackedPileup that is free of mapping quality zero reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no MQ0 reads in the pileup.
     *
     * @return
     */
    public ReadBackedPileup getPileupWithoutMappingQualityZeroReads() {
        SampleSplitReadBackedPileup filteredPileup = new SampleSplitReadBackedPileup(loc);

        for (Map.Entry<String,ReadBackedPileup> entry: pileupBySample.entrySet()) {
            String sampleName = entry.getKey();
            ReadBackedPileup pileup = entry.getValue();
            filteredPileup.addPileup(sampleName,pileup.getPileupWithoutMappingQualityZeroReads());
        }

        return filteredPileup;
    }

    /**
     * Returns subset of this pileup that contains only bases with quality >= minBaseQ, coming from
     * reads with mapping qualities >= minMapQ. This method allocates and returns a new instance of ReadBackedPileup.
     * @param minBaseQ
     * @param minMapQ
     * @return
     */
    @Override
    public ReadBackedPileup getBaseAndMappingFilteredPileup( int minBaseQ, int minMapQ ) {
        SampleSplitReadBackedPileup filteredPileup = new SampleSplitReadBackedPileup(loc);

        for (Map.Entry<String,ReadBackedPileup> entry: pileupBySample.entrySet()) {
            String sampleName = entry.getKey();
            ReadBackedPileup pileup = entry.getValue();
            filteredPileup.addPileup(sampleName,pileup.getBaseAndMappingFilteredPileup(minBaseQ,minMapQ));
        }

        return filteredPileup;
    }

    /** Returns subset of this pileup that contains only bases with quality >= minBaseQ.
     * This method allocates and returns a new instance of ReadBackedPileup.
     * @param minBaseQ
     * @return
     */
    @Override
    public ReadBackedPileup getBaseFilteredPileup( int minBaseQ ) {
        return getBaseAndMappingFilteredPileup(minBaseQ, -1);
    }

    /** Returns subset of this pileup that contains only bases coming from reads with mapping quality >= minMapQ.
     * This method allocates and returns a new instance of ReadBackedPileup.
     * @param minMapQ
     * @return
     */
    @Override
    public ReadBackedPileup getMappingFilteredPileup( int minMapQ ) {
        return getBaseAndMappingFilteredPileup(-1, minMapQ);
    }

    /**
     * Returns a pileup randomly downsampled to the desiredCoverage.
     *
     * @param desiredCoverage
     * @return
     */
    @Override
    public ReadBackedPileup getDownsampledPileup(int desiredCoverage) {
        if ( size() <= desiredCoverage )
            return this;

        // randomly choose numbers corresponding to positions in the reads list
        Random generator = new Random();
        TreeSet<Integer> positions = new TreeSet<Integer>();
        for ( int i = 0; i < desiredCoverage; /* no update */ ) {
            if ( positions.add(generator.nextInt(size)) )
                i++;
        }

        SampleSplitReadBackedPileup downsampledPileup = new SampleSplitReadBackedPileup(getLocation());
        int current = 0;

        for(Map.Entry<String,ReadBackedPileup> entry: this.pileupBySample.entrySet()) {
            String sampleName = entry.getKey();
            ReadBackedPileup pileup = entry.getValue();

            List<PileupElement> filteredPileup = new ArrayList<PileupElement>();
            for(PileupElement p: pileup) {
                if(positions.contains(current))
                    filteredPileup.add(p);
            }

            if(!filteredPileup.isEmpty())
                downsampledPileup.addPileup(sampleName,new UnifiedReadBackedPileup(loc,filteredPileup));

            current++;
        }

        return downsampledPileup;
    }

    @Override
    public ReadBackedPileup getPileupForSample(String sampleName) {
        return pileupBySample.containsKey(sampleName) ? pileupBySample.get(sampleName) : null;
    }

    // --------------------------------------------------------
    //
    // iterators
    //
    // --------------------------------------------------------

    /**
     * The best way to access PileupElements where you only care about the bases and quals in the pileup.
     *
     * for (PileupElement p : this) { doSomething(p); }
     *
     * Provides efficient iteration of the data.
     *
     * @return
     */
    @Override
    public Iterator<PileupElement> iterator() {
        return new MergingPileupElementIterator(pileupBySample);
    }

    // todo -- why is this public?
    @Override
    public IterableIterator<ExtendedPileupElement> extendedForeachIterator() {
        ArrayList<ExtendedPileupElement> x = new ArrayList<ExtendedPileupElement>(size());
        int i = 0;
        Iterator<PileupElement> iterator = new MergingPileupElementIterator(pileupBySample);
        while(iterator.hasNext()) {
            PileupElement pile = iterator.next();
            x.add(new ExtendedPileupElement(pile.getRead(), pile.getOffset(), i++, this));
        }

        return new IterableIterator<ExtendedPileupElement>(x.iterator());
    }

    /**
     * Simple useful routine to count the number of deletion bases in this pileup
     *
     * @return
     */
    @Override
    public int getNumberOfDeletions() {
        return nDeletions;
    }

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

    /**
     * Get counts of A, C, G, T in order, which returns a int[4] vector with counts according
     * to BaseUtils.simpleBaseToBaseIndex for each base.
     *
     * @return
     */
    @Override
    public int[] getBaseCounts() {
        int[] counts = new int[4];
        for(ReadBackedPileup pileup: pileupBySample.values()) {
            int[] countsBySample = pileup.getBaseCounts();
            for(int i = 0; i < counts.length; i++)
                counts[i] += countsBySample[i];
        }
        return counts;        
    }

    /**
     * Somewhat expensive routine that returns true if any base in the pileup has secondary bases annotated
     * @return
     */
    @Override
    public boolean hasSecondaryBases() {
        boolean hasSecondaryBases = false;
        for(ReadBackedPileup pileup: pileupBySample.values())
            hasSecondaryBases |= pileup.hasSecondaryBases();
        return hasSecondaryBases;
    }

    @Override
    public String getPileupString(char ref) {
        // In the pileup format, each line represents a genomic position, consisting of chromosome name,
        // coordinate, reference base, read bases, read qualities and alignment mapping qualities.
        return String.format("%s %s %c %s %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                ref,                                                     // reference base
                new String(getBases()),
                getQualsString());
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
        for ( PileupElement pile : this ) { reads.add(pile.getRead()); }
        return reads;
    }

    /**
     * Returns a list of the offsets in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    @Override
    public List<Integer> getOffsets() {
        List<Integer> offsets = new ArrayList<Integer>(size());
        for ( PileupElement pile : this ) { offsets.add(pile.getOffset()); }
        return offsets;
    }

    /**
     * Returns an array of the bases in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    @Override
    public byte[] getBases() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getBase(); }
        return v;
    }

    /**
     * Returns an array of the secondary bases in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    @Override
    public byte[] getSecondaryBases() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getSecondBase(); }
        return v;
    }

     /**
     * Returns an array of the quals in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
     @Override
     public byte[] getQuals() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getQual(); }
        return v;
     }

    /**
     * Get an array of the mapping qualities
     * @return
     */
    @Override
    public byte[] getMappingQuals() {
       byte[] v = new byte[size()];
       for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = (byte)pile.getRead().getMappingQuality(); }
       return v;
    }

    private String getQualsString() {
        return UnifiedReadBackedPileup.quals2String(getQuals());
    }
}
