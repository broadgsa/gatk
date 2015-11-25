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

package org.broadinstitute.gatk.utils.pileup;

import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.fragments.FragmentCollection;
import org.broadinstitute.gatk.utils.fragments.FragmentUtils;
import org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.BaseUtils;

import java.util.*;

public class ReadBackedPileupImpl implements ReadBackedPileup {
    protected final GenomeLoc loc;
    protected final PileupElementTracker<PileupElement> pileupElementTracker;

    private final static int UNINITIALIZED_CACHED_INT_VALUE = -1;

    /**
     * Different then number of elements due to reduced reads
     */
    private int depthOfCoverage = UNINITIALIZED_CACHED_INT_VALUE;
    private int nDeletions = UNINITIALIZED_CACHED_INT_VALUE;            // cached value of the number of deletions
    private int nMQ0Reads = UNINITIALIZED_CACHED_INT_VALUE;             // cached value of the number of MQ0 reads

    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.  This pileup will contain a list, in order of the reads, of the piled bases at
     * reads[i] for all i in offsets.  Does not make a copy of the data, so it's not safe to
     * go changing the reads.
     *
     * @param loc     The genome loc to associate reads wotj
     * @param reads
     * @param offsets
     */
    public ReadBackedPileupImpl(GenomeLoc loc, List<GATKSAMRecord> reads, List<Integer> offsets) {
        this.loc = loc;
        this.pileupElementTracker = readsOffsets2Pileup(reads, offsets);
    }


    /**
     * Create a new version of a read backed pileup at loc without any aligned reads
     */
    public ReadBackedPileupImpl(GenomeLoc loc) {
        this(loc, new UnifiedPileupElementTracker<PileupElement>());
    }

    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.  This lower level constructure assumes pileup is well-formed and merely keeps a
     * pointer to pileup.  Don't go changing the data in pileup.
     */
    public ReadBackedPileupImpl(GenomeLoc loc, List<PileupElement> pileup) {
        if (loc == null) throw new ReviewedGATKException("Illegal null genomeloc in ReadBackedPileup");
        if (pileup == null) throw new ReviewedGATKException("Illegal null pileup in ReadBackedPileup");

        this.loc = loc;
        this.pileupElementTracker = new UnifiedPileupElementTracker<PileupElement>(pileup);
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     *
     * @param loc
     * @param pileup
     */
    @Deprecated
    public ReadBackedPileupImpl(GenomeLoc loc, List<PileupElement> pileup, int size, int nDeletions, int nMQ0Reads) {
        this(loc, pileup);
    }

    protected ReadBackedPileupImpl(GenomeLoc loc, PileupElementTracker<PileupElement> tracker) {
        this.loc = loc;
        this.pileupElementTracker = tracker;
    }

    public ReadBackedPileupImpl(GenomeLoc loc, Map<String, ReadBackedPileupImpl> pileupsBySample) {
        this.loc = loc;
        PerSamplePileupElementTracker<PileupElement> tracker = new PerSamplePileupElementTracker<PileupElement>();
        for (Map.Entry<String, ReadBackedPileupImpl> pileupEntry : pileupsBySample.entrySet()) {
            tracker.addElements(pileupEntry.getKey(), pileupEntry.getValue().pileupElementTracker);
        }
        this.pileupElementTracker = tracker;
    }

    public ReadBackedPileupImpl(GenomeLoc loc, List<GATKSAMRecord> reads, int offset) {
        this.loc = loc;
        this.pileupElementTracker = readsOffsets2Pileup(reads, offset);
    }

    /**
     * Helper routine for converting reads and offset lists to a PileupElement list.
     *
     * @param reads
     * @param offsets
     * @return
     */
    private PileupElementTracker<PileupElement> readsOffsets2Pileup(List<GATKSAMRecord> reads, List<Integer> offsets) {
        if (reads == null) throw new ReviewedGATKException("Illegal null read list in UnifiedReadBackedPileup");
        if (offsets == null) throw new ReviewedGATKException("Illegal null offsets list in UnifiedReadBackedPileup");
        if (reads.size() != offsets.size())
            throw new ReviewedGATKException("Reads and offset lists have different sizes!");

        UnifiedPileupElementTracker<PileupElement> pileup = new UnifiedPileupElementTracker<PileupElement>();
        for (int i = 0; i < reads.size(); i++) {
            GATKSAMRecord read = reads.get(i);
            int offset = offsets.get(i);
            pileup.add(createNewPileupElement(read, offset)); // only used to create fake pileups for testing so ancillary information is not important
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
    private PileupElementTracker<PileupElement> readsOffsets2Pileup(List<GATKSAMRecord> reads, int offset) {
        if (reads == null) throw new ReviewedGATKException("Illegal null read list in UnifiedReadBackedPileup");
        if (offset < 0) throw new ReviewedGATKException("Illegal offset < 0 UnifiedReadBackedPileup");

        UnifiedPileupElementTracker<PileupElement> pileup = new UnifiedPileupElementTracker<PileupElement>();
        for (GATKSAMRecord read : reads) {
            pileup.add(createNewPileupElement(read, offset)); // only used to create fake pileups for testing so ancillary information is not important
        }

        return pileup;
    }

    protected ReadBackedPileupImpl createNewPileup(GenomeLoc loc, PileupElementTracker<PileupElement> tracker) {
        return new ReadBackedPileupImpl(loc, tracker);
    }

    protected PileupElement createNewPileupElement(GATKSAMRecord read, int offset) {
        return LocusIteratorByState.createPileupForReadAndOffset(read, offset);
    }    
    
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
    public ReadBackedPileupImpl getPileupWithoutDeletions() {
        if (getNumberOfDeletions() > 0) {
            if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
                PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
                PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

                for (final String sample : tracker.getSamples()) {
                    PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                    ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getPileupWithoutDeletions();
                    filteredTracker.addElements(sample, pileup.pileupElementTracker);
                }
                return createNewPileup(loc, filteredTracker);

            } else {
                UnifiedPileupElementTracker<PileupElement> tracker = (UnifiedPileupElementTracker<PileupElement>) pileupElementTracker;
                UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();

                for (PileupElement p : tracker) {
                    if (!p.isDeletion()) {
                        filteredTracker.add(p);
                    }
                }
                return createNewPileup(loc, filteredTracker);
            }
        } else {
            return this;
        }
    }

    /**
     * Returns a new ReadBackedPileup where only one read from an overlapping read
     * pair is retained.  If the two reads in question disagree to their basecall,
     * neither read is retained.  If they agree on the base, the read with the higher
     * base quality observation is retained
     *
     * @return the newly filtered pileup
     */
    @Override
    public ReadBackedPileup getOverlappingFragmentFilteredPileup() {
        return getOverlappingFragmentFilteredPileup(true, true);
    }

    /**
     * Returns a new ReadBackedPileup where only one read from an overlapping read
     * pair is retained.  If discardDiscordant and the two reads in question disagree to their basecall,
     * neither read is retained.  Otherwise, the read with the higher
     * quality (base or mapping, depending on baseQualNotMapQual) observation is retained
     *
     * @return the newly filtered pileup
     */
    @Override
    public ReadBackedPileupImpl getOverlappingFragmentFilteredPileup(boolean discardDiscordant, boolean baseQualNotMapQual) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getOverlappingFragmentFilteredPileup(discardDiscordant, baseQualNotMapQual);
                filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }
            return createNewPileup(loc, filteredTracker);
        } else {
            Map<String, PileupElement> filteredPileup = new HashMap<String, PileupElement>();

            for (PileupElement p : pileupElementTracker) {
                String readName = p.getRead().getReadName();

                // if we've never seen this read before, life is good
                if (!filteredPileup.containsKey(readName)) {
                    filteredPileup.put(readName, p);
                } else {
                    PileupElement existing = filteredPileup.get(readName);

                    // if the reads disagree at this position, throw them both out.  Otherwise
                    // keep the element with the higher quality score
                    if (discardDiscordant && existing.getBase() != p.getBase()) {
                        filteredPileup.remove(readName);
                    } else {
                        if (baseQualNotMapQual) {
                            if (existing.getQual() < p.getQual())
                                filteredPileup.put(readName, p);
                        }
                        else {
                            if (existing.getMappingQual() < p.getMappingQual())
                                filteredPileup.put(readName, p);
                        }
                    }
                }
            }

            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();
            for (PileupElement filteredElement : filteredPileup.values())
                filteredTracker.add(filteredElement);

            return createNewPileup(loc, filteredTracker);
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
    public ReadBackedPileupImpl getPileupWithoutMappingQualityZeroReads() {
        if (getNumberOfMappingQualityZeroReads() > 0) {
            if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
                PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
                PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

                for (final String sample : tracker.getSamples()) {
                    PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                    ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getPileupWithoutMappingQualityZeroReads();
                    filteredTracker.addElements(sample, pileup.pileupElementTracker);
                }
                return createNewPileup(loc, filteredTracker);

            } else {
                UnifiedPileupElementTracker<PileupElement> tracker = (UnifiedPileupElementTracker<PileupElement>) pileupElementTracker;
                UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();

                for (PileupElement p : tracker) {
                    if (p.getRead().getMappingQuality() > 0) {
                        filteredTracker.add(p);
                    }
                }
                return createNewPileup(loc, filteredTracker);
            }
        } else {
            return this;
        }
    }

    public ReadBackedPileupImpl getPositiveStrandPileup() {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getPositiveStrandPileup();
                filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }
            return createNewPileup(loc, filteredTracker);
        } else {
            UnifiedPileupElementTracker<PileupElement> tracker = (UnifiedPileupElementTracker<PileupElement>) pileupElementTracker;
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();

            for (PileupElement p : tracker) {
                if (!p.getRead().getReadNegativeStrandFlag()) {
                    filteredTracker.add(p);
                }
            }
            return createNewPileup(loc, filteredTracker);
        }
    }

    /**
     * Gets the pileup consisting of only reads on the negative strand.
     *
     * @return A read-backed pileup consisting only of reads on the negative strand.
     */
    public ReadBackedPileupImpl getNegativeStrandPileup() {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getNegativeStrandPileup();
                filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }
            return createNewPileup(loc, filteredTracker);
        } else {
            UnifiedPileupElementTracker<PileupElement> tracker = (UnifiedPileupElementTracker<PileupElement>) pileupElementTracker;
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();

            for (PileupElement p : tracker) {
                if (p.getRead().getReadNegativeStrandFlag()) {
                    filteredTracker.add(p);
                }
            }
            return createNewPileup(loc, filteredTracker);
        }
    }

    /**
     * Gets a pileup consisting of all those elements passed by a given filter.
     *
     * @param filter Filter to use when testing for elements.
     * @return a pileup without the given filtered elements.
     */
    public ReadBackedPileupImpl getFilteredPileup(PileupElementFilter filter) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getFilteredPileup(filter);
                filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }

            return createNewPileup(loc, filteredTracker);
        } else {
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();

            for (PileupElement p : pileupElementTracker) {
                if (filter.allow(p))
                    filteredTracker.add(p);
            }

            return createNewPileup(loc, filteredTracker);
        }
    }

    /**
     * Returns subset of this pileup that contains only bases with quality >= minBaseQ, coming from
     * reads with mapping qualities >= minMapQ. This method allocates and returns a new instance of ReadBackedPileup.
     *
     * @param minBaseQ
     * @param minMapQ
     * @return
     */
    @Override
    public ReadBackedPileupImpl getBaseAndMappingFilteredPileup(int minBaseQ, int minMapQ) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getBaseAndMappingFilteredPileup(minBaseQ, minMapQ);
                filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }

            return createNewPileup(loc, filteredTracker);
        } else {
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();

            for (PileupElement p : pileupElementTracker) {
                if (p.getRead().getMappingQuality() >= minMapQ && (p.isDeletion() || p.getQual() >= minBaseQ)) {
                    filteredTracker.add(p);
                }
            }

            return createNewPileup(loc, filteredTracker);
        }
    }

    /**
     * Returns subset of this pileup that contains only bases with quality >= minBaseQ.
     * This method allocates and returns a new instance of ReadBackedPileup.
     *
     * @param minBaseQ
     * @return
     */
    @Override
    public ReadBackedPileup getBaseFilteredPileup(int minBaseQ) {
        return getBaseAndMappingFilteredPileup(minBaseQ, -1);
    }

    /**
     * Returns subset of this pileup that contains only bases coming from reads with mapping quality >= minMapQ.
     * This method allocates and returns a new instance of ReadBackedPileup.
     *
     * @param minMapQ
     * @return
     */
    @Override
    public ReadBackedPileup getMappingFilteredPileup(int minMapQ) {
        return getBaseAndMappingFilteredPileup(-1, minMapQ);
    }

    /**
     * Gets a list of the read groups represented in this pileup.
     *
     * @return
     */
    @Override
    public Collection<String> getReadGroups() {
        Set<String> readGroups = new HashSet<String>();
        for (PileupElement pileupElement : this)
            readGroups.add(pileupElement.getRead().getReadGroup().getReadGroupId());
        return readGroups;
    }

    /**
     * Gets the pileup for a given read group.  Horrendously inefficient at this point.
     *
     * @param targetReadGroupId Identifier for the read group.
     * @return A read-backed pileup containing only the reads in the given read group.
     */
    @Override
    public ReadBackedPileupImpl getPileupForReadGroup(String targetReadGroupId) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getPileupForReadGroup(targetReadGroupId);
                if (pileup != null)
                    filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        } else {
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();
            for (PileupElement p : pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if (targetReadGroupId != null) {
                    if (read.getReadGroup() != null && targetReadGroupId.equals(read.getReadGroup().getReadGroupId()))
                        filteredTracker.add(p);
                } else {
                    if (read.getReadGroup() == null || read.getReadGroup().getReadGroupId() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        }
    }

    /**
     * Gets the pileup for a set of read groups.  Horrendously inefficient at this point.
     *
     * @param rgSet List of identifiers for the read groups.
     * @return A read-backed pileup containing only the reads in the given read groups.
     */
    @Override
    public ReadBackedPileupImpl getPileupForReadGroups(final HashSet<String> rgSet) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getPileupForReadGroups(rgSet);
                if (pileup != null)
                    filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        } else {
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();
            for (PileupElement p : pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if (rgSet != null && !rgSet.isEmpty()) {
                    if (read.getReadGroup() != null && rgSet.contains(read.getReadGroup().getReadGroupId()))
                        filteredTracker.add(p);
                } else {
                    if (read.getReadGroup() == null || read.getReadGroup().getReadGroupId() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        }
    }

    @Override
    public ReadBackedPileupImpl getPileupForLane(String laneID) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                ReadBackedPileupImpl pileup = createNewPileup(loc, perSampleElements).getPileupForLane(laneID);
                if (pileup != null)
                    filteredTracker.addElements(sample, pileup.pileupElementTracker);
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        } else {
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();
            for (PileupElement p : pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if (laneID != null) {
                    if (read.getReadGroup() != null &&
                            (read.getReadGroup().getReadGroupId().startsWith(laneID + ".")) ||   // lane is the same, but sample identifier is different
                            (read.getReadGroup().getReadGroupId().equals(laneID)))               // in case there is no sample identifier, they have to be exactly the same
                        filteredTracker.add(p);
                } else {
                    if (read.getReadGroup() == null || read.getReadGroup().getReadGroupId() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        }
    }

    public Collection<String> getSamples() {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            return new HashSet<String>(tracker.getSamples());
        } else {
            Collection<String> sampleNames = new HashSet<String>();
            for (PileupElement p : this) {
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
     * TODO: delete this once the experimental downsampler stabilizes
     *
     * @param desiredCoverage
     * @return
     */
    @Override
    public ReadBackedPileup getDownsampledPileup(int desiredCoverage) {
        if (getNumberOfElements() <= desiredCoverage)
            return this;

        // randomly choose numbers corresponding to positions in the reads list
        TreeSet<Integer> positions = new TreeSet<Integer>();
        for (int i = 0; i < desiredCoverage; /* no update */) {
            if (positions.add(Utils.getRandomGenerator().nextInt(getNumberOfElements())))
                i++;
        }

        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PerSamplePileupElementTracker<PileupElement> filteredTracker = new PerSamplePileupElementTracker<PileupElement>();


            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);

                int current = 0;
                UnifiedPileupElementTracker<PileupElement> filteredPileup = new UnifiedPileupElementTracker<PileupElement>();
                for (PileupElement p : perSampleElements) {
                    if (positions.contains(current))
                        filteredPileup.add(p);
                    current++;

                }
                filteredTracker.addElements(sample, filteredPileup);
            }

            return createNewPileup(loc, filteredTracker);
        } else {
            UnifiedPileupElementTracker<PileupElement> tracker = (UnifiedPileupElementTracker<PileupElement>) pileupElementTracker;
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();

            Iterator positionIter = positions.iterator();

            while (positionIter.hasNext()) {
                int nextReadToKeep = (Integer) positionIter.next();
                filteredTracker.add(tracker.get(nextReadToKeep));
            }

            return createNewPileup(getLocation(), filteredTracker);
        }
    }

    @Override
    public ReadBackedPileup getPileupForSamples(Collection<String> sampleNames) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PileupElementTracker<PileupElement> filteredElements = tracker.getElements(sampleNames);
            return filteredElements != null ? createNewPileup(loc, filteredElements) : null;
        } else {
            HashSet<String> hashSampleNames = new HashSet<String>(sampleNames);                                         // to speed up the "contains" access in the for loop
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();
            for (PileupElement p : pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if (sampleNames != null) {                                                                              // still checking on sampleNames because hashSampleNames will never be null. And empty means something else.
                    if (read.getReadGroup() != null && hashSampleNames.contains(read.getReadGroup().getSample()))
                        filteredTracker.add(p);
                } else {
                    if (read.getReadGroup() == null || read.getReadGroup().getSample() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        }
    }

    @Override
    public Map<String, ReadBackedPileup> getPileupsForSamples(Collection<String> sampleNames) {
        Map<String, ReadBackedPileup> result = new HashMap<String, ReadBackedPileup>();
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            for (String sample : sampleNames) {
                PileupElementTracker<PileupElement> filteredElements = tracker.getElements(sample);
                if (filteredElements != null)
                    result.put(sample, createNewPileup(loc, filteredElements));
            }
        } else {
            Map<String, UnifiedPileupElementTracker<PileupElement>> trackerMap = new HashMap<String, UnifiedPileupElementTracker<PileupElement>>();

            for (String sample : sampleNames) {                                                                         // initialize pileups for each sample
                UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();
                trackerMap.put(sample, filteredTracker);
            }
            for (PileupElement p : pileupElementTracker) {                                                                         // go through all pileup elements only once and add them to the respective sample's pileup
                GATKSAMRecord read = p.getRead();
                if (read.getReadGroup() != null) {
                    String sample = read.getReadGroup().getSample();
                    UnifiedPileupElementTracker<PileupElement> tracker = trackerMap.get(sample);
                    if (tracker != null)                                                                                // we only add the pileup the requested samples. Completely ignore the rest
                        tracker.add(p);
                }
            }
            for (Map.Entry<String, UnifiedPileupElementTracker<PileupElement>> entry : trackerMap.entrySet())                      // create the ReadBackedPileup for each sample
                result.put(entry.getKey(), createNewPileup(loc, entry.getValue()));
        }
        return result;
    }


    @Override
    public ReadBackedPileup getPileupForSample(String sampleName) {
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            PileupElementTracker<PileupElement> filteredElements = tracker.getElements(sampleName);
            return filteredElements != null ? createNewPileup(loc, filteredElements) : null;
        } else {
            UnifiedPileupElementTracker<PileupElement> filteredTracker = new UnifiedPileupElementTracker<PileupElement>();
            for (PileupElement p : pileupElementTracker) {
                GATKSAMRecord read = p.getRead();
                if (sampleName != null) {
                    if (read.getReadGroup() != null && sampleName.equals(read.getReadGroup().getSample()))
                        filteredTracker.add(p);
                } else {
                    if (read.getReadGroup() == null || read.getReadGroup().getSample() == null)
                        filteredTracker.add(p);
                }
            }
            return filteredTracker.size() > 0 ? createNewPileup(loc, filteredTracker) : null;
        }
    }

    // --------------------------------------------------------
    //
    // iterators
    //
    // --------------------------------------------------------

    /**
     * The best way to access PileupElements where you only care about the bases and quals in the pileup.
     * <p/>
     * for (PileupElement p : this) { doSomething(p); }
     * <p/>
     * Provides efficient iteration of the data.
     *
     * @return
     */
    @Override
    public Iterator<PileupElement> iterator() {
        return new Iterator<PileupElement>() {
            private final Iterator<PileupElement> wrappedIterator = pileupElementTracker.iterator();

            public boolean hasNext() {
                return wrappedIterator.hasNext();
            }

            public PileupElement next() {
                return wrappedIterator.next();
            }

            public void remove() {
                throw new UnsupportedOperationException("Cannot remove from a pileup element iterator");
            }
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
        if ( nDeletions == UNINITIALIZED_CACHED_INT_VALUE ) {
            nDeletions = 0;
            for (PileupElement p : pileupElementTracker.unorderedIterable() ) {
                if (p.isDeletion()) {
                    nDeletions++;
                }
            }
        }
        return nDeletions;
    }

    @Override
    public int getNumberOfMappingQualityZeroReads() {
        if ( nMQ0Reads == UNINITIALIZED_CACHED_INT_VALUE ) {
            nMQ0Reads = 0;

            for (PileupElement p : pileupElementTracker.unorderedIterable()) {
                if (p.getRead().getMappingQuality() == 0) {
                    nMQ0Reads++;
                }
            }
        }

        return nMQ0Reads;
    }

    /**
     * @return the number of physical elements in this pileup
     */
    @Override
    public int getNumberOfElements() {
        return pileupElementTracker.size();
    }

    /**
     * @return the number of abstract elements in this pileup
     */
    @Override
    public int depthOfCoverage() {
        if (depthOfCoverage == UNINITIALIZED_CACHED_INT_VALUE) {
            depthOfCoverage = pileupElementTracker.size();
        }
        return depthOfCoverage;
    }

    /**
     * @return true if there are 0 elements in the pileup, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return getNumberOfElements() == 0;
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

        // TODO -- can be optimized with .unorderedIterable()
        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;
            for (final String sample : tracker.getSamples()) {
                int[] countsBySample = createNewPileup(loc, tracker.getElements(sample)).getBaseCounts();
                for (int i = 0; i < counts.length; i++)
                    counts[i] += countsBySample[i];
            }
        } else {
            for (PileupElement pile : this) {
                // skip deletion sites
                if (!pile.isDeletion()) {
                    int index = BaseUtils.simpleBaseToBaseIndex((char) pile.getBase());
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
     *
     * @return
     */
    @Override
    public List<GATKSAMRecord> getReads() {
        List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>(getNumberOfElements());
        for (PileupElement pile : this) {
            reads.add(pile.getRead());
        }
        return reads;
    }

    @Override
    public int getNumberOfDeletionsAfterThisElement() {
        int count = 0;
        for (PileupElement p : pileupElementTracker.unorderedIterable()) {
            if (p.isBeforeDeletionStart())
                count++;
        }
        return count;
    }

    @Override
    public int getNumberOfInsertionsAfterThisElement() {
        int count = 0;
        for (PileupElement p : pileupElementTracker.unorderedIterable()) {
            if (p.isBeforeInsertion())
                count++;
        }
        return count;

    }
    /**
     * Returns a list of the offsets in this pileup. Note this call costs O(n) and allocates fresh lists each time
     *
     * @return
     */
    @Override
    public List<Integer> getOffsets() {
        List<Integer> offsets = new ArrayList<Integer>(getNumberOfElements());
        for (PileupElement pile : pileupElementTracker.unorderedIterable()) {
            offsets.add(pile.getOffset());
        }
        return offsets;
    }

    /**
     * Returns an array of the bases in this pileup. Note this call costs O(n) and allocates fresh array each time
     *
     * @return
     */
    @Override
    public byte[] getBases() {
        byte[] v = new byte[getNumberOfElements()];
        int pos = 0;
        for (PileupElement pile : pileupElementTracker) {
            v[pos++] = pile.getBase();
        }
        return v;
    }

    /**
     * Returns an array of the quals in this pileup. Note this call costs O(n) and allocates fresh array each time
     *
     * @return
     */
    @Override
    public byte[] getQuals() {
        byte[] v = new byte[getNumberOfElements()];
        int pos = 0;
        for (PileupElement pile : pileupElementTracker) {
            v[pos++] = pile.getQual();
        }
        return v;
    }

    /**
     * Get an array of the mapping qualities
     *
     * @return
     */
    @Override
    public int[] getMappingQuals() {
        final int[] v = new int[getNumberOfElements()];
        int pos = 0;
        for ( final PileupElement pile : pileupElementTracker ) {
            v[pos++] = pile.getRead().getMappingQuality();
        }
        return v;
    }

    static String quals2String(byte[] quals) {
        StringBuilder qualStr = new StringBuilder();
        for (int qual : quals) {
            qual = Math.min(qual, 63);              // todo: fixme, this isn't a good idea
            char qualChar = (char) (33 + qual);     // todo: warning, this is illegal for qual > 63
            qualStr.append(qualChar);
        }

        return qualStr.toString();
    }

    private String getQualsString() {
        return quals2String(getQuals());
    }

    /**
     * Returns a new ReadBackedPileup that is sorted by start coordinate of the reads.
     *
     * @return
     */
    @Override
    public ReadBackedPileup getStartSortedPileup() {

        final TreeSet<PileupElement> sortedElements = new TreeSet<PileupElement>(new Comparator<PileupElement>() {
            @Override
            public int compare(PileupElement element1, PileupElement element2) {
                final int difference = element1.getRead().getAlignmentStart() - element2.getRead().getAlignmentStart();
                return difference != 0 ? difference : element1.getRead().getReadName().compareTo(element2.getRead().getReadName());
            }
        });

        if (pileupElementTracker instanceof PerSamplePileupElementTracker) {
            PerSamplePileupElementTracker<PileupElement> tracker = (PerSamplePileupElementTracker<PileupElement>) pileupElementTracker;

            for (final String sample : tracker.getSamples()) {
                PileupElementTracker<PileupElement> perSampleElements = tracker.getElements(sample);
                for (PileupElement pile : perSampleElements)
                    sortedElements.add(pile);
            }
        }
        else {
            UnifiedPileupElementTracker<PileupElement> tracker = (UnifiedPileupElementTracker<PileupElement>) pileupElementTracker;
            for (PileupElement pile : tracker)
                sortedElements.add(pile);
        }

        UnifiedPileupElementTracker<PileupElement> sortedTracker = new UnifiedPileupElementTracker<PileupElement>();
        for (PileupElement pile : sortedElements)
            sortedTracker.add(pile);

        return createNewPileup(loc, sortedTracker);
    }

    @Override
    public FragmentCollection<PileupElement> toFragments() {
        return FragmentUtils.create(this);
    }

    @Override
    public ReadBackedPileup copy() {
        return new ReadBackedPileupImpl(loc, pileupElementTracker.copy());
    }
}
