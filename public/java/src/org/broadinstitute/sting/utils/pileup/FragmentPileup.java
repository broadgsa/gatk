package org.broadinstitute.sting.utils.pileup;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

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
 * User: depristo
 * Date: 3/26/11
 * Time: 10:09 PM
 */
public class FragmentPileup {
    final Collection<PileupElement> oneReadPile;
    final Collection<TwoReadPileupElement> twoReadPile = new ArrayList<TwoReadPileupElement>();

    /**
     * Create a new Fragment-based pileup from the standard read-based pileup
     * @param pileup
     */
    public FragmentPileup(ReadBackedPileup pileup) {
        Map<String, PileupElement> nameMap = new HashMap<String, PileupElement>();

        // build an initial map, grabbing all of the multi-read fragments
        for ( PileupElement p : pileup ) {
            String readName = p.getRead().getReadName();

            PileupElement pe1 = nameMap.get(readName);
            if ( pe1 != null ) {
                // assumes we have at most 2 reads per fragment
                twoReadPile.add(new TwoReadPileupElement(pe1, p));
                nameMap.remove(readName);
            } else {
                nameMap.put(readName, p);
            }
        }

        // now set the one Read pile to the values in the nameMap with only a single read
        oneReadPile = nameMap.values();
    }

    /**
     * Gets the pileup elements containing two reads, in no particular order
     *
     * @return
     */
    public Collection<TwoReadPileupElement> getTwoReadPileup() {
        return twoReadPile;
    }

    /**
     * Gets the pileup elements containing one read, in no particular order
     *
     * @return
     */
    public Collection<PileupElement> getOneReadPileup() {
        return oneReadPile;
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
