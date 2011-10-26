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

    protected enum FragmentMatchingAlgorithm {
        ORIGINAL,
        skipNonOverlapping,
    }

    /**
     * Create a new Fragment-based pileup from the standard read-based pileup
     * @param pileup
     */
    public FragmentPileup(ReadBackedPileup pileup) {
        skipNonOverlapping(pileup);
    }

    /** For performance testing only */
    protected FragmentPileup(ReadBackedPileup pileup, FragmentMatchingAlgorithm algorithm) {
        switch ( algorithm ) {
            case ORIGINAL: oldSlowCalculation(pileup); break;
            case skipNonOverlapping: skipNonOverlapping(pileup); break;
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

    private final void skipNonOverlapping(final ReadBackedPileup pileup) {
        Map<String, PileupElement> nameMap = null;

        // build an initial map, grabbing all of the multi-read fragments
        for ( final PileupElement p : pileup ) {
            final SAMRecord read = p.getRead();
            final int mateStart = read.getMateAlignmentStart();

            if ( mateStart == 0 || mateStart > read.getAlignmentEnd() ) {
                // if we know that this read won't overlap its mate, or doesn't have one, jump out early
                if ( oneReadPile == null ) oneReadPile = new ArrayList<PileupElement>(pileup.size()); // lazy init
                oneReadPile.add(p);
            } else {
                // the read might overlap it's mate, or is the rightmost read of a pair
                final String readName = p.getRead().getReadName();
                final PileupElement pe1 = nameMap == null ? null : nameMap.get(readName);
                if ( pe1 != null ) {
                    // assumes we have at most 2 reads per fragment
                    if ( twoReadPile == null ) twoReadPile = new ArrayList<TwoReadPileupElement>(); // lazy init
                    twoReadPile.add(new TwoReadPileupElement(pe1, p));
                    nameMap.remove(readName);
                } else {
                    if ( nameMap == null ) nameMap = new HashMap<String, PileupElement>(pileup.size()); // lazy init
                    nameMap.put(readName, p);
                }
            }
        }

        // add all of the reads that are potentially overlapping but whose mate never showed
        // up to the oneReadPile
        if ( nameMap != null && ! nameMap.isEmpty() ) {
            if ( oneReadPile == null )
                oneReadPile = nameMap.values();
            else
                oneReadPile.addAll(nameMap.values());
        }
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
