package org.broadinstitute.sting.utils.pileup;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.*;

/**
 * An easy to access fragment-based pileup, which contains two separate pileups.  The first
 * is a regular collection of PileupElements containing all of the reads in the original RBP
 * that uniquely info about a fragment.  The second are TwoReadPileupElements that, as the
 * name suggests, contain two reads that are sequenced from the same underlying fragment.
 *
 * Based on the original code by E. Banks
 *
 * TODO -- technically we could generalize this code to support a pseudo-duplicate marking
 * TODO -- algorithm that could collect all duplicates into a single super pileup element
 *
 * Oct 21: note that the order of the oneReadPileup and twoReadPileups are not
 * defined.  The algorithms that produce these lists are in fact producing
 * lists of Pileup elements *NOT* sorted by alignment start position of the underlying
 * reads.
 *
 * User: depristo
 * Date: 3/26/11
 * Time: 10:09 PM
 */
public class FragmentPileup {
    Collection<PileupElement> oneReadPile = null;
    Collection<TwoReadPileupElement> twoReadPile = null;

    public enum FragmentMatchingAlgorithm {
        ORIGINAL,
        FAST_V1,
        skipNonOverlapping,
        skipNonOverlappingNotLazy
    }

    /**
     * Create a new Fragment-based pileup from the standard read-based pileup
     * @param pileup
     */
    public FragmentPileup(ReadBackedPileup pileup) {
        //oldSlowCalculation(pileup);
        fastNewCalculation(pileup);
    }

    protected FragmentPileup(ReadBackedPileup pileup, FragmentMatchingAlgorithm algorithm) {
        switch ( algorithm ) {
            case ORIGINAL: oldSlowCalculation(pileup); break;
            case FAST_V1: fastNewCalculation(pileup); break;
            case skipNonOverlapping: skipNonOverlapping(pileup); break;
            case skipNonOverlappingNotLazy: skipNonOverlappingNotLazy(pileup); break;
        }
    }

    private final void oldSlowCalculation(final ReadBackedPileup pileup) {
        final Map<String, PileupElement> nameMap = new HashMap<String, PileupElement>(pileup.size());

        // build an initial map, grabbing all of the multi-read fragments
        for ( final PileupElement p : pileup ) {
            final String readName = p.getRead().getReadName();

            final PileupElement pe1 = nameMap.get(readName);
            if ( pe1 != null ) {
                // assumes we have at most 2 reads per fragment
                if ( twoReadPile == null ) twoReadPile = new ArrayList<TwoReadPileupElement>();
                twoReadPile.add(new TwoReadPileupElement(pe1, p));
                nameMap.remove(readName);
            } else {
                nameMap.put(readName, p);
            }
        }

        oneReadPile = nameMap.values();
    }

    /**
     * @param pileup
     */
    private final void fastNewCalculation(final ReadBackedPileup pileup) {
        Map<String, PileupElement> nameMap = null; // lazy initialization

        for ( final PileupElement p : pileup ) {
            final SAMRecord read = p.getRead();

            switch (ReadUtils.readMightOverlapMate(read) ) {
                // we know for certain this read doesn't have an overlapping mate
                case NO: {
                    addToOnePile(p);
                    break;
                }

                // we know that we overlap our mate, so put the read in the nameMap in
                // case our mate shows up
                case LEFT_YES: {
                    nameMap = addToNameMap(nameMap, p);
                    break;
                }

                // read starts at the same position, so we are looking at either the first or
                // the second read.  In the first, add it to the map, in the second grab it
                // from the map and create a fragment
                case SAME_START: {
                    final PileupElement pe1 = getFromNameMap(nameMap, p);
                    if ( pe1 != null ) {
                        addToTwoPile(pe1, p);
                        nameMap.remove(p.getRead().getReadName());
                    } else {
                        nameMap = addToNameMap(nameMap, p);
                    }
                    break;
                }

                // in this case we need to see if our mate is already present, and if so
                // grab the read from the list
                case RIGHT_MAYBE: {
                    final PileupElement pe1 = getFromNameMap(nameMap, p);
                    if ( pe1 != null ) {
                        addToTwoPile(pe1, p);
                        nameMap.remove(p.getRead().getReadName());
                    } else {
                        addToOnePile(p);
                    }
                    break;
                }
            }
        }

        if ( nameMap != null && ! nameMap.isEmpty() ) {
            if ( oneReadPile == null )
                oneReadPile = nameMap.values();
            else
                oneReadPile.addAll(nameMap.values());
        }
    }

    /**
     * @param pileup
     */
    private final void skipNonOverlappingNotLazy(final ReadBackedPileup pileup) {
        oneReadPile = new ArrayList<PileupElement>(pileup.size());
        twoReadPile = new ArrayList<TwoReadPileupElement>();
        final Map<String, PileupElement> nameMap = new HashMap<String, PileupElement>(pileup.size());

        // build an initial map, grabbing all of the multi-read fragments
        for ( final PileupElement p : pileup ) {
            // if we know that this read won't overlap its mate, or doesn't have one, jump out early
            final SAMRecord read = p.getRead();
            final int mateStart = read.getMateAlignmentStart();
            if ( mateStart == 0 || mateStart > read.getAlignmentEnd() ) {
                oneReadPile.add(p);
            } else {
                final String readName = p.getRead().getReadName();
                final PileupElement pe1 = nameMap.get(readName);
                if ( pe1 != null ) {
                    // assumes we have at most 2 reads per fragment
                    twoReadPile.add(new TwoReadPileupElement(pe1, p));
                    nameMap.remove(readName);
                } else {
                    nameMap.put(readName, p);
                }
            }
        }

        oneReadPile.addAll(nameMap.values());
    }

    private final void skipNonOverlapping(final ReadBackedPileup pileup) {
        Map<String, PileupElement> nameMap = null;

        // build an initial map, grabbing all of the multi-read fragments
        for ( final PileupElement p : pileup ) {
            // if we know that this read won't overlap its mate, or doesn't have one, jump out early
            final SAMRecord read = p.getRead();
            final int mateStart = read.getMateAlignmentStart();
            if ( mateStart == 0 || mateStart > read.getAlignmentEnd() ) {
                if ( oneReadPile == null ) oneReadPile = new ArrayList<PileupElement>(pileup.size());
                oneReadPile.add(p);
            } else {
                final String readName = p.getRead().getReadName();
                final PileupElement pe1 = nameMap == null ? null : nameMap.get(readName);
                if ( pe1 != null ) {
                    // assumes we have at most 2 reads per fragment
                    if ( twoReadPile == null ) twoReadPile = new ArrayList<TwoReadPileupElement>();
                    twoReadPile.add(new TwoReadPileupElement(pe1, p));
                    nameMap.remove(readName);
                } else {
                    nameMap = addToNameMap(nameMap, p);
                }
            }
        }

        if ( oneReadPile == null )
            oneReadPile = nameMap == null ? Collections.<PileupElement>emptyList() : nameMap.values();
        else if ( nameMap != null )
            oneReadPile.addAll(nameMap.values());
    }

    private final Map<String, PileupElement> addToNameMap(Map<String, PileupElement> map, final PileupElement p) {
        if ( map == null ) map = new HashMap<String, PileupElement>();
        map.put(p.getRead().getReadName(), p);
        return map;
    }

    private final PileupElement getFromNameMap(Map<String, PileupElement> map, final PileupElement p) {
        return map == null ? null : map.get(p.getRead().getReadName());
    }


    private final void addToOnePile(final PileupElement p) {
        if ( oneReadPile == null ) oneReadPile = new ArrayList<PileupElement>();
        oneReadPile.add(p);
    }

    private final void addToTwoPile(final PileupElement p1, final PileupElement p2) {
        // assumes we have at most 2 reads per fragment
        if ( twoReadPile == null ) twoReadPile = new ArrayList<TwoReadPileupElement>();
        twoReadPile.add(new TwoReadPileupElement(p1, p2));
    }

    /**
     * Gets the pileup elements containing two reads, in no particular order
     *
     * @return
     */
    public Collection<TwoReadPileupElement> getTwoReadPileup() {
        return twoReadPile == null ? Collections.<TwoReadPileupElement>emptyList() : twoReadPile;
    }

    /**
     * Gets the pileup elements containing one read, in no particular order
     *
     * @return
     */
    public Collection<PileupElement> getOneReadPileup() {
        return oneReadPile == null ? Collections.<PileupElement>emptyList() : oneReadPile;
    }

    /**
     * Useful helper class to represent a full read pair at a position
     *
     * User: ebanks, depristo
     * Date: Jan 10, 2011
     */
    public static class TwoReadPileupElement {
        final protected PileupElement PE1, PE2;

        /**
         * Creates a fragment element that contains both ends of a paired end read
         * @param PE1
         * @param PE2
         */
        public TwoReadPileupElement(PileupElement PE1, PileupElement PE2) {
            this.PE1 = PE1;
            this.PE2 = PE2;
        }

        /** Returns the first pileup element -- never null */
        public PileupElement getFirst() { return PE1; }

        /** Returns the second read in this fragment element.  May be null */
        public PileupElement getSecond() { return PE2; }
    }
}
