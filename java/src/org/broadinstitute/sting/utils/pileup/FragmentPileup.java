package org.broadinstitute.sting.utils.pileup;

import java.util.*;

/**
 * An easy to access fragment-based pileup.  new FragmentPileup(RBP) creates one, and you
 * can either iterate over or get the collection of PerFragmentPileupElements.
 *
 * Based on the original code by E. Banks
 *
 * User: depristo
 * Date: 3/26/11
 * Time: 10:09 PM
 */
public class FragmentPileup implements Iterable<FragmentPileup.PerFragmentPileupElement> {
    final Collection<PerFragmentPileupElement> fragments = new ArrayList<PerFragmentPileupElement>();

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
                fragments.add(new PerFragmentPileupElement(pe1, p));
                nameMap.remove(readName);
            } else {
                nameMap.put(readName, p);
            }
        }

        // now go through the values in the nameMap to get the fragments with only a single read
        for ( PileupElement p : nameMap.values() )
            fragments.add(new PerFragmentPileupElement(p));
    }

    /**
     * Gets the fragments, in no particular order
     *
     * @return
     */
    public Collection<PerFragmentPileupElement> getFragments() {
        return fragments;
    }

    /**
     * Returns an iterator over the fragments.  No specific order of fragments is assumed
     * @return
     */
    public Iterator<PerFragmentPileupElement> iterator() {
        return fragments.iterator();
    }

    /**
     * Useful helper class to represent a full read pair at a position
     *
     * User: ebanks
     * Date: Jan 10, 2011
     */
    public static class PerFragmentPileupElement {
        protected PileupElement PE1 = null, PE2 = null;

        /**
         * Creates a fragment element that only contains a single read
         * @param PE
         */
        public PerFragmentPileupElement(PileupElement PE) {
            PE1 = PE;
        }

        /**
         * Creates a fragment element that contains both ends of a paired end read
         * @param PE1
         * @param PE2
         */
        public PerFragmentPileupElement(PileupElement PE1, PileupElement PE2) {
            this.PE1 = PE1;
            this.PE2 = PE2;
        }

        /** Returns the first pileup element -- never null */
        public PileupElement getFirst() { return PE1; }

        /** Is there a second read in this fragment element? */
        public boolean hasSecond() { return PE2 != null; }

        /** Returns the second read in this fragment element.  May be null */
        public PileupElement getSecond() { return PE2; }
    }
}
