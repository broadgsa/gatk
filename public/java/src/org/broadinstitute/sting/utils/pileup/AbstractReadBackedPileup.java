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

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.fragments.FragmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * A generic implementation of read-backed pileups.
 *
 * @author mhanna
 * @version 0.1
 */
public abstract class AbstractReadBackedPileup<RBP extends AbstractReadBackedPileup<RBP,PE>,PE extends PileupElement> implements ReadBackedPileup {
    protected final GenomeLoc loc;
    protected final PileupElementTracker<PE> pileupElementTracker;

    protected int size = 0;                   // cached value of the size of the pileup
    protected int abstractSize = -1;           // cached value of the abstract size of the pileup
    protected int nDeletions = 0;             // cached value of the number of deletions
    protected int nMQ0Reads = 0;              // cached value of the number of MQ0 reads

    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.  This pileup will contain a list, in order of the reads, of the piled bases at
     * reads[i] for all i in offsets.  Does not make a copy of the data, so it's not safe to
     * go changing the reads.
     *
     * @param loc
     * @param reads
     * @param offsets
     */
    public AbstractReadBackedPileup(GenomeLoc loc, List<GATKSAMRecord> reads, List<Integer> offsets ) {
        this.loc = loc;
        this.pileupElementTracker = readsOffsets2Pileup(reads,offsets);
    }

    public AbstractReadBackedPileup(GenomeLoc loc, List<GATKSAMRecord> reads, int offset ) {
        this.loc = loc;
        this.pileupElementTracker = readsOffsets2Pileup(reads,offset);
    }

    /**
     * Create a new version of a read backed pileup at loc without any aligned reads
     *
     */
    public AbstractReadBackedPileup(GenomeLoc loc) {
        this(loc, new UnifiedPileupElementTracker<PE>());
    }

    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.  This lower level constructure assumes pileup is well-formed and merely keeps a
     * pointer to pileup.  Don't go changing the data in pileup.
     *
     */
    public AbstractReadBackedPileup(GenomeLoc loc, List<PE> pileup) {
        if ( loc == null ) throw new ReviewedStingException("Illegal null genomeloc in ReadBackedPileup");
        if ( pileup == null ) throw new ReviewedStingException("Illegal null pileup in ReadBackedPileup");

        this.loc = loc;
        this.pileupElementTracker = new UnifiedPileupElementTracker<PE>(pileup);
        calculateCachedData();
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc
     * @param pileup
     */
    public AbstractReadBackedPileup(GenomeLoc loc, List<PE> pileup, int size, int nDeletions, int nMQ0Reads) {
        if ( loc == null ) throw new ReviewedStingException("Illegal null genomeloc in UnifiedReadBackedPileup");
        if ( pileup == null ) throw new ReviewedStingException("Illegal null pileup in UnifiedReadBackedPileup");

        this.loc = loc;
        this.pileupElementTracker = new UnifiedPileupElementTracker<PE>(pileup);
        this.size = size;
        this.nDeletions = nDeletions;
        this.nMQ0Reads = nMQ0Reads;
    }


    protected AbstractReadBackedPileup(GenomeLoc loc, PileupElementTracker<PE> tracker) {
        this.loc = loc;
        this.pileupElementTracker = tracker;
        calculateCachedData();
    }

    protected AbstractReadBackedPileup(GenomeLoc loc, Map<String,? extends AbstractReadBackedPileup<RBP,PE>> pileupsBySample) {
        this.loc = loc;
        PerSamplePileupElementTracker<PE> tracker = new PerSamplePileupElementTracker<PE>();
        for(Map.Entry<String,? extends AbstractReadBackedPileup<RBP,PE>> pileupEntry: pileupsBySample.entrySet()) {
            tracker.addElements(pileupEntry.getKey(),pileupEntry.getValue().pileupElementTracker);
            addPileupToCumulativeStats(pileupEntry.getValue());
        }
        this.pileupElementTracker = tracker;
    }

    /**
     * Calculate cached sizes, nDeletion, and base counts for the pileup.  This calculation is done upfront,
     * so you pay the cost at the start, but it's more efficient to do this rather than pay the cost of calling
     * sizes, nDeletion, etc. over and over potentially.
     */
    protected void calculateCachedData() {
        size = 0;
        nDeletions = 0;
        nMQ0Reads = 0;

        for ( PileupElement p : pileupElementTracker ) {
            size++;
            if ( p.isDeletion() ) {
                nDeletions++;
            }
            if ( p.getRead().getMappingQuality() == 0 ) {
                nMQ0Reads++;
            }
        }
    }

    protected void calculateAbstractSize() {
        abstractSize = 0;
        for ( PileupElement p : pileupElementTracker ) {
            abstractSize += p.getRepresentativeCount();
        }
    }

    protected void addPileupToCumulativeStats(AbstractReadBackedPileup<RBP,PE> pileup) {
        size += pileup.getNumberOfElements();
        abstractSize += pileup.depthOfCoverage();
        nDeletions += pileup.getNumberOfDeletions();
        nMQ0Reads += pileup.getNumberOfMappingQualityZeroReads();
    }

    /**
     * Helper routine for converting reads and offset lists to a PileupElement list.
     *
     * @param reads
     * @param offsets
     * @return
     */
    private PileupElementTracker<PE> readsOffsets2Pileup(List<GATKSAMRecord> reads, List<Integer> offsets ) {
        if ( reads == null ) throw new ReviewedStingException("Illegal null read list in UnifiedReadBackedPileup");
        if ( offsets == null ) throw new ReviewedStingException("Illegal null offsets list in UnifiedReadBackedPileup");
        if ( reads.size() != offsets.size() ) throw new ReviewedStingException("Reads and offset lists have different sizes!");

        UnifiedPileupElementTracker<PE> pileup = new UnifiedPileupElementTracker<PE>();
        for ( int i = 0; i < reads.size(); i++ ) {
            pileup.add(createNewPileupElement(reads.get(i),offsets.get(i)));
        }

        return pileup;
    }

    /**
     * Helper routine for converting reads and a single offset to a PileupElement list.
     *
     * @param reads
     * @param offset
     * @return
     */
    private PileupElementTracker<PE> readsOffsets2Pileup(List<GATKSAMRecord> reads, int offset ) {
        if ( reads == null ) throw new ReviewedStingException("Illegal null read list in UnifiedReadBackedPileup");
        if ( offset < 0 ) throw new ReviewedStingException("Illegal offset < 0 UnifiedReadBackedPileup");

        UnifiedPileupElementTracker<PE> pileup = new UnifiedPileupElementTracker<PE>();
        for ( int i = 0; i < reads.size(); i++ ) {
            pileup.add(createNewPileupElement( reads.get(i), offset ));
        }

        return pileup;
    }

    protected abstract AbstractReadBackedPileup<RBP,PE> createNewPileup(GenomeLoc loc, PileupElementTracker<PE> pileupElementTracker);
    protected abstract PE createNewPileupElement(GATKSAMRecord read, int offset);

    // --------------------------------------------------------
    //
    // Special 'constructors'
    //
    // --------------------------------------------------------

    /**
     * Returns a new ReadBackedPileup that is free of deletion spanning reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no deletions in the pileup.
     *
     * @return
     */
    @Override
    public RBP getPileupWithoutDeletions() {
        if ( getNumberOfDeletions() > 0 ) {
            if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
                PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
                PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

                for(final String sample: tracker.getSamples()) {
                    PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                    AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getPileupWithoutDeletions();
                    filteredTracker.addElements(sample,pileup.pileupElementTracker);
                }
                return (RBP)createNewPileup(loc,filteredTracker);

            }
            else {
                UnifiedPileupElementTracker<PE> tracker = (UnifiedPileupElementTracker<PE>)pileupElementTracker;
                UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();

                for ( PE p : tracker ) {
                    if ( !p.isDeletion() ) {
                        filteredTracker.add(p);
                    }
                }
                return (RBP)createNewPileup(loc, filteredTracker);
            }
        } else {
            return (RBP)this;
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
    public RBP getOverlappingFragmentFilteredPileup() {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getOverlappingFragmentFilteredPileup();
                filteredTracker.addElements(sample,pileup.pileupElementTracker);
            }
            return (RBP)createNewPileup(loc,filteredTracker);
        }
        else {
            Map<String,PE> filteredPileup = new HashMap<String, PE>();

            for ( PE p : pileupElementTracker ) {
                String readName = p.getRead().getReadName();

                // if we've never seen this read before, life is good
                if (!filteredPileup.containsKey(readName)) {
                    filteredPileup.put(readName, p);
                } else {
                    PileupElement existing = filteredPileup.get(readName);

                    // if the reads disagree at this position, throw them both out.  Otherwise
                    // keep the element with the higher quality score
                    if (existing.getBase() != p.getBase()) {
                        filteredPileup.remove(readName);
                    } else {
                        if (existing.getQual() < p.getQual()) {
                            filteredPileup.put(readName, p);
                        }
                    }
                }
            }

            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();
            for(PE filteredElement: filteredPileup.values())
                filteredTracker.add(filteredElement);

            return (RBP)createNewPileup(loc,filteredTracker);
        }
    }


    /**
     * Returns a new ReadBackedPileup that is free of mapping quality zero reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no MQ0 reads in the pileup.
     *
     * @return
     */
    @Override
    public RBP getPileupWithoutMappingQualityZeroReads() {
        if ( getNumberOfMappingQualityZeroReads() > 0 ) {
            if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
                PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
                PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

                for(final String sample: tracker.getSamples()) {
                    PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                    AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getPileupWithoutMappingQualityZeroReads();
                    filteredTracker.addElements(sample,pileup.pileupElementTracker);
                }
                return (RBP)createNewPileup(loc,filteredTracker);

            }
            else {
                UnifiedPileupElementTracker<PE> tracker = (UnifiedPileupElementTracker<PE>)pileupElementTracker;
                UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();

                for ( PE p : tracker ) {
                    if ( p.getRead().getMappingQuality() > 0 ) {
                        filteredTracker.add(p);
                    }
                }
                return (RBP)createNewPileup(loc, filteredTracker);
            }
        } else {
            return (RBP)this;
        }
    }

    public RBP getPositiveStrandPileup() {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getPositiveStrandPileup();
                filteredTracker.addElements(sample,pileup.pileupElementTracker);
            }
            return (RBP)createNewPileup(loc,filteredTracker);
        }
        else {
            UnifiedPileupElementTracker<PE> tracker = (UnifiedPileupElementTracker<PE>)pileupElementTracker;
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();

            for ( PE p : tracker ) {
                if ( !p.getRead().getReadNegativeStrandFlag() ) {
                    filteredTracker.add(p);
                }
            }
            return (RBP)createNewPileup(loc, filteredTracker);
        }
    }

    /**
     * Gets the pileup consisting of only reads on the negative strand.
     * @return A read-backed pileup consisting only of reads on the negative strand.
     */
    public RBP getNegativeStrandPileup() {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getNegativeStrandPileup();
                filteredTracker.addElements(sample,pileup.pileupElementTracker);
            }
            return (RBP)createNewPileup(loc,filteredTracker);
        }
        else {
            UnifiedPileupElementTracker<PE> tracker = (UnifiedPileupElementTracker<PE>)pileupElementTracker;
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();

            for ( PE p : tracker ) {
                if ( p.getRead().getReadNegativeStrandFlag() ) {
                    filteredTracker.add(p);
                }
            }
            return (RBP)createNewPileup(loc, filteredTracker);
        }
    }

    /**
     * Gets a pileup consisting of all those elements passed by a given filter.
     * @param filter Filter to use when testing for elements.
     * @return a pileup without the given filtered elements.
     */
    public RBP getFilteredPileup(PileupElementFilter filter) {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getFilteredPileup(filter);
                filteredTracker.addElements(sample,pileup.pileupElementTracker);
            }

            return (RBP)createNewPileup(loc,filteredTracker);
        }
        else {
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();

            for ( PE p : pileupElementTracker ) {
                if( filter.allow(p) )
                    filteredTracker.add(p);
            }

            return (RBP)createNewPileup(loc, filteredTracker);
        }
    }

    /** Returns subset of this pileup that contains only bases with quality >= minBaseQ, coming from
     * reads with mapping qualities >= minMapQ. This method allocates and returns a new instance of ReadBackedPileup.
     * @param minBaseQ
     * @param minMapQ
     * @return
     */
    @Override
    public RBP getBaseAndMappingFilteredPileup( int minBaseQ, int minMapQ ) {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getBaseAndMappingFilteredPileup(minBaseQ,minMapQ);
                filteredTracker.addElements(sample,pileup.pileupElementTracker);
            }

            return (RBP)createNewPileup(loc,filteredTracker);
        }
        else {
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();

            for ( PE p : pileupElementTracker ) {
                if ( p.getRead().getMappingQuality() >= minMapQ &&
                        (p.isDeletion() ||
                                ((p instanceof ExtendedEventPileupElement) && ((ExtendedEventPileupElement)p).getType() == ExtendedEventPileupElement.Type.NOEVENT) ||
                                p.getQual() >= minBaseQ) ) {
                    filteredTracker.add(p);
                }
            }

            return (RBP)createNewPileup(loc, filteredTracker);
        }
    }

    /** Returns subset of this pileup that contains only bases with quality >= minBaseQ.
     * This method allocates and returns a new instance of ReadBackedPileup.
     * @param minBaseQ
     * @return
     */
    @Override
    public RBP getBaseFilteredPileup( int minBaseQ ) {
        return getBaseAndMappingFilteredPileup(minBaseQ, -1);
    }

    /** Returns subset of this pileup that contains only bases coming from reads with mapping quality >= minMapQ.
     * This method allocates and returns a new instance of ReadBackedPileup.
     * @param minMapQ
     * @return
     */
    @Override
    public RBP getMappingFilteredPileup( int minMapQ ) {
        return getBaseAndMappingFilteredPileup(-1, minMapQ);
    }

    /**
     * Gets a list of the read groups represented in this pileup.
     * @return
     */
    @Override
    public Collection<String> getReadGroups() {
        Set<String> readGroups = new HashSet<String>();
        for(PileupElement pileupElement: this)
            readGroups.add(pileupElement.getRead().getReadGroup().getReadGroupId());
        return readGroups;
    }

    /**
     * Gets the pileup for a given read group.  Horrendously inefficient at this point.
     * @param targetReadGroupId Identifier for the read group.
     * @return A read-backed pileup containing only the reads in the given read group.
     */
    @Override
    public RBP getPileupForReadGroup(String targetReadGroupId) {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getPileupForReadGroup(targetReadGroupId);
                if(pileup != null)
                    filteredTracker.addElements(sample,pileup.pileupElementTracker);
            }
            return filteredTracker.size()>0 ? (RBP)createNewPileup(loc,filteredTracker) : null;
        }
        else {
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();
            for(PE p: pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if(targetReadGroupId != null) {
                    if(read.getReadGroup() != null && targetReadGroupId.equals(read.getReadGroup().getReadGroupId()))
                        filteredTracker.add(p);
                }
                else {
                    if(read.getReadGroup() == null || read.getReadGroup().getReadGroupId() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size()>0 ? (RBP)createNewPileup(loc,filteredTracker) : null;
        }
    }

    @Override
    public RBP getPileupForLane(String laneID) {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);
                AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements).getPileupForLane(laneID);
                if(pileup != null)
                    filteredTracker.addElements(sample,pileup.pileupElementTracker);
            }
            return filteredTracker.size()>0 ? (RBP)createNewPileup(loc,filteredTracker) : null;
        }
        else {
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();
            for(PE p: pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if(laneID != null) {
                    if(read.getReadGroup() != null &&
                       (read.getReadGroup().getReadGroupId().startsWith(laneID + ".")) ||   // lane is the same, but sample identifier is different
                       (read.getReadGroup().getReadGroupId().equals(laneID)))               // in case there is no sample identifier, they have to be exactly the same
                        filteredTracker.add(p);
                }
                else {
                    if(read.getReadGroup() == null || read.getReadGroup().getReadGroupId() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size()>0 ? (RBP)createNewPileup(loc,filteredTracker) : null;
        }
    }

    public Collection<String> getSamples() {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            return new HashSet<String>(tracker.getSamples());
        }
        else {
            Collection<String> sampleNames = new HashSet<String>();
            for(PileupElement p: this) {
                GATKSAMRecord read = p.getRead();
                String sampleName = read.getReadGroup() != null ? read.getReadGroup().getSample() : null;
                sampleNames.add(sampleName);
            }
            return sampleNames;
        }
    }

    /**
     * Returns a pileup randomly downsampled to the desiredCoverage.
     *
     * @param desiredCoverage
     * @return
     */
    @Override
    public RBP getDownsampledPileup(int desiredCoverage) {
        if ( getNumberOfElements() <= desiredCoverage )
            return (RBP)this;

        // randomly choose numbers corresponding to positions in the reads list
        TreeSet<Integer> positions = new TreeSet<Integer>();
        for ( int i = 0; i < desiredCoverage; /* no update */ ) {
            if ( positions.add(GenomeAnalysisEngine.getRandomGenerator().nextInt(size)) )
                i++;
        }

        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PerSamplePileupElementTracker<PE> filteredTracker = new PerSamplePileupElementTracker<PE>();

            int current = 0;

            for(final String sample: tracker.getSamples()) {
                PileupElementTracker<PE> perSampleElements = tracker.getElements(sample);

                List<PileupElement> filteredPileup = new ArrayList<PileupElement>();
                for(PileupElement p: perSampleElements) {
                    if(positions.contains(current))
                        filteredPileup.add(p);
                }

                if(!filteredPileup.isEmpty()) {
                    AbstractReadBackedPileup<RBP,PE> pileup = createNewPileup(loc,perSampleElements);
                    filteredTracker.addElements(sample,pileup.pileupElementTracker);
                }

                current++;
            }

            return (RBP)createNewPileup(loc,filteredTracker);
        }
        else {
            UnifiedPileupElementTracker<PE> tracker = (UnifiedPileupElementTracker<PE>)pileupElementTracker;
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();

            Iterator positionIter = positions.iterator();

            while ( positionIter.hasNext() ) {
                int nextReadToKeep = (Integer)positionIter.next();
                filteredTracker.add(tracker.get(nextReadToKeep));
            }

            return (RBP)createNewPileup(getLocation(), filteredTracker);
        }
    }

    @Override
    public RBP getPileupForSamples(Collection<String> sampleNames) {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PileupElementTracker<PE> filteredElements = tracker.getElements(sampleNames);
            return filteredElements != null ? (RBP)createNewPileup(loc,filteredElements) : null;
        }
        else {
            HashSet<String> hashSampleNames = new HashSet<String>(sampleNames);    // to speed up the "contains" access in the for loop
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();
            for(PE p: pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if(sampleNames != null) {                                          // still checking on sampleNames because hashSampleNames will never be null. And empty means something else.
                    if(read.getReadGroup() != null && hashSampleNames.contains(read.getReadGroup().getSample()))
                        filteredTracker.add(p);
                }
                else {
                    if(read.getReadGroup() == null || read.getReadGroup().getSample() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size()>0 ? (RBP)createNewPileup(loc,filteredTracker) : null;
        }
    }


    @Override
    public RBP getPileupForSample(String sampleName) {
        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            PileupElementTracker<PE> filteredElements = tracker.getElements(sampleName);
            return filteredElements != null ? (RBP)createNewPileup(loc,filteredElements) : null;
        }
        else {
            UnifiedPileupElementTracker<PE> filteredTracker = new UnifiedPileupElementTracker<PE>();
            for(PE p: pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if(sampleName != null) {
                    if(read.getReadGroup() != null && sampleName.equals(read.getReadGroup().getSample()))
                        filteredTracker.add(p);
                }
                else {
                    if(read.getReadGroup() == null || read.getReadGroup().getSample() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size()>0 ? (RBP)createNewPileup(loc,filteredTracker) : null;
        }
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
        return new Iterator<PileupElement>() {
            private final Iterator<PE> wrappedIterator = pileupElementTracker.iterator();

            public boolean hasNext() { return wrappedIterator.hasNext(); }
            public PileupElement next() { return wrappedIterator.next(); }
            public void remove() { throw new UnsupportedOperationException("Cannot remove from a pileup element iterator"); }
        };
    }

    /**
     * The best way to access PileupElements where you only care not only about bases and quals in the pileup
     * but also need access to the index of the pileup element in the pile.
     *
     * for (ExtendedPileupElement p : this) { doSomething(p); }
     *
     * Provides efficient iteration of the data.
     *
     * @return
     */

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
     * @return the number of physical elements in this pileup
     */
    @Override
    public int getNumberOfElements() {
        return size;
    }

    /**
     * @return the number of abstract elements in this pileup
     */
    @Override
    public int depthOfCoverage() {
        if ( abstractSize == -1 )
            calculateAbstractSize();
        return abstractSize;
    }

    /**
     * @return true if there are 0 elements in the pileup, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return size==0;
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

        if(pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PE> tracker = (PerSamplePileupElementTracker<PE>)pileupElementTracker;
            for(final String sample: tracker.getSamples()) {
                int[] countsBySample = createNewPileup(loc,tracker.getElements(sample)).getBaseCounts();
                for(int i = 0; i < counts.length; i++)
                    counts[i] += countsBySample[i];
            }
        }
        else {
            for ( PileupElement pile : this ) {
                // skip deletion sites
                if ( ! pile.isDeletion() ) {
                    int index = BaseUtils.simpleBaseToBaseIndex((char)pile.getBase());
                    if (index != -1)
                        counts[index]++;
                }
            }
        }

        return counts;
    }

    @Override
    public String getPileupString(Character ref) {
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
    public List<GATKSAMRecord> getReads() {
        List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>(getNumberOfElements());
        for ( PileupElement pile : this ) { reads.add(pile.getRead()); }
        return reads;
    }

    /**
     * Returns a list of the offsets in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    @Override
    public List<Integer> getOffsets() {
        List<Integer> offsets = new ArrayList<Integer>(getNumberOfElements());
        for ( PileupElement pile : this ) { offsets.add(pile.getOffset()); }
        return offsets;
    }

    /**
     * Returns an array of the bases in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    @Override
    public byte[] getBases() {
        byte[] v = new byte[getNumberOfElements()];
        int pos = 0;
        for ( PileupElement pile : pileupElementTracker ) { v[pos++] = pile.getBase(); }
        return v;
    }

    /**
     * Returns an array of the quals in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    @Override
    public byte[] getQuals() {
        byte[] v = new byte[getNumberOfElements()];
        int pos = 0;
        for ( PileupElement pile : pileupElementTracker ) { v[pos++] = pile.getQual(); }
        return v;
    }

    /**
     * Get an array of the mapping qualities
     * @return
     */
    @Override
    public byte[] getMappingQuals() {
        byte[] v = new byte[getNumberOfElements()];
        int pos = 0;
        for ( PileupElement pile : pileupElementTracker ) { v[pos++] = (byte)pile.getRead().getMappingQuality(); }
        return v;
    }

    static String quals2String( byte[] quals ) {
        StringBuilder qualStr = new StringBuilder();
        for ( int qual : quals ) {
            qual = Math.min(qual, 63);              // todo: fixme, this isn't a good idea
            char qualChar = (char) (33 + qual);     // todo: warning, this is illegal for qual > 63
            qualStr.append(qualChar);
        }

        return qualStr.toString();
    }

    private String getQualsString() {
        return quals2String(getQuals());
    }

    @Override
    public FragmentCollection<PileupElement> toFragments() {
        return FragmentUtils.create(this);
    }
}


