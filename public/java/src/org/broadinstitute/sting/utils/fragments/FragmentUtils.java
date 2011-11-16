package org.broadinstitute.sting.utils.fragments;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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
public class FragmentUtils {
    private FragmentUtils() {} // private constructor

    /**
     * A getter function that takes an Object of type T and returns its associated SAMRecord.
     *
     * Allows us to write a generic T -> Fragment algorithm that works with any object containing
     * a read.
     *
     * @param <T>
     */
    public interface ReadGetter<T> {
        public GATKSAMRecord get(T object);
    }

    /** Identify getter for SAMRecords themselves */
    private final static ReadGetter<GATKSAMRecord> SamRecordGetter = new ReadGetter<GATKSAMRecord>() {
        @Override public GATKSAMRecord get(final GATKSAMRecord object) { return object; }
    };

    /** Gets the SAMRecord in a PileupElement */
    private final static ReadGetter<PileupElement> PileupElementGetter = new ReadGetter<PileupElement>() {
        @Override public GATKSAMRecord get(final PileupElement object) { return object.getRead(); }
    };


    /**
     * Generic algorithm that takes an iterable over T objects, a getter routine to extract the reads in T,
     * and returns a FragmentCollection that contains the T objects whose underlying reads either overlap (or
     * not) with their mate pairs.
     *
     * @param readContainingObjects
     * @param nElements
     * @param getter
     * @param <T>
     * @return
     */
    private final static <T> FragmentCollection<T> create(Iterable<T> readContainingObjects, int nElements, ReadGetter<T> getter) {
        Collection<T> singletons = null;
        Collection<List<T>> overlapping = null;
        Map<String, T> nameMap = null;

        int lastStart = -1;

        // build an initial map, grabbing all of the multi-read fragments
        for ( final T p : readContainingObjects ) {
            final SAMRecord read = getter.get(p);

            if ( read.getAlignmentStart() < lastStart ) {
                throw new IllegalArgumentException(String.format(
                        "FragmentUtils.create assumes that the incoming objects are ordered by " +
                                "SAMRecord alignment start, but saw a read %s with alignment start " +
                                "%d before the previous start %d", read.getSAMString(), read.getAlignmentStart(), lastStart));
            }
            lastStart = read.getAlignmentStart();

            final int mateStart = read.getMateAlignmentStart();
            if ( mateStart == 0 || mateStart > read.getAlignmentEnd() ) {
                // if we know that this read won't overlap its mate, or doesn't have one, jump out early
                if ( singletons == null ) singletons = new ArrayList<T>(nElements); // lazy init
                singletons.add(p);
            } else {
                // the read might overlap it's mate, or is the rightmost read of a pair
                final String readName = read.getReadName();
                final T pe1 = nameMap == null ? null : nameMap.get(readName);
                if ( pe1 != null ) {
                    // assumes we have at most 2 reads per fragment
                    if ( overlapping == null ) overlapping = new ArrayList<List<T>>(); // lazy init
                    overlapping.add(Arrays.asList(pe1, p));
                    nameMap.remove(readName);
                } else {
                    if ( nameMap == null ) nameMap = new HashMap<String, T>(nElements); // lazy init
                    nameMap.put(readName, p);
                }
            }
        }

        // add all of the reads that are potentially overlapping but whose mate never showed
        // up to the oneReadPile
        if ( nameMap != null && ! nameMap.isEmpty() ) {
            if ( singletons == null )
                singletons = nameMap.values();
            else
                singletons.addAll(nameMap.values());
        }

        return new FragmentCollection<T>(singletons, overlapping);
    }

    public final static FragmentCollection<PileupElement> create(ReadBackedPileup rbp) {
        return create(rbp, rbp.getNumberOfElements(), PileupElementGetter);
    }

    public final static FragmentCollection<GATKSAMRecord> create(List<GATKSAMRecord> reads) {
        return create(reads, reads.size(), SamRecordGetter);
    }

}
