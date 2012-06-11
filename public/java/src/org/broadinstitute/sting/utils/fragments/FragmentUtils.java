package org.broadinstitute.sting.utils.fragments;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.walkers.bqsr.EventType;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
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

    public final static List<GATKSAMRecord> mergeOverlappingPairedFragments( List<GATKSAMRecord> overlappingPair ) {
        final byte MIN_QUAL_BAD_OVERLAP = 16;
        if( overlappingPair.size() != 2 ) { throw new ReviewedStingException("Found overlapping pair with " + overlappingPair.size() + " reads, but expecting exactly 2."); }

        GATKSAMRecord firstRead = overlappingPair.get(0);
        GATKSAMRecord secondRead = overlappingPair.get(1);
        if( !(secondRead.getUnclippedStart() <= firstRead.getUnclippedEnd() && secondRead.getUnclippedStart() >= firstRead.getUnclippedStart() && secondRead.getUnclippedEnd() >= firstRead.getUnclippedEnd()) ) {
            firstRead = overlappingPair.get(1); // swap them
            secondRead = overlappingPair.get(0);
        }
        if( !(secondRead.getUnclippedStart() <= firstRead.getUnclippedEnd() && secondRead.getUnclippedStart() >= firstRead.getUnclippedStart() && secondRead.getUnclippedEnd() >= firstRead.getUnclippedEnd()) ) {
            return overlappingPair; // can't merge them, yet:  AAAAAAAAAAA-BBBBBBBBBBB-AAAAAAAAAAAAAA, B is contained entirely inside A
        }
        if( firstRead.getCigarString().contains("I") || firstRead.getCigarString().contains("D") || secondRead.getCigarString().contains("I") || secondRead.getCigarString().contains("D") ) {
            return overlappingPair; // fragments contain indels so don't merge them
        }

        final Pair<Integer, Boolean> pair = ReadUtils.getReadCoordinateForReferenceCoordinate(firstRead, secondRead.getSoftStart());

        final int firstReadStop = ( pair.getSecond() ? pair.getFirst() + 1 : pair.getFirst() );
        final int numBases = firstReadStop + secondRead.getReadLength();
        final byte[] bases = new byte[numBases];
        final byte[] quals = new byte[numBases];
        final byte[] insertionQuals = new byte[numBases];
        final byte[] deletionQuals = new byte[numBases];
        final byte[] firstReadBases = firstRead.getReadBases();
        final byte[] firstReadQuals = firstRead.getBaseQualities();
        final byte[] secondReadBases = secondRead.getReadBases();
        final byte[] secondReadQuals = secondRead.getBaseQualities();

        for(int iii = 0; iii < firstReadStop; iii++) {
            bases[iii] = firstReadBases[iii];
            quals[iii] = firstReadQuals[iii];
        }
        for(int iii = firstReadStop; iii < firstRead.getReadLength(); iii++) {
            if( firstReadQuals[iii] > MIN_QUAL_BAD_OVERLAP && secondReadQuals[iii-firstReadStop] > MIN_QUAL_BAD_OVERLAP && firstReadBases[iii] != secondReadBases[iii-firstReadStop] ) {
                return overlappingPair;// high qual bases don't match exactly, probably indel in only one of the fragments, so don't merge them
            }
            if( firstReadQuals[iii] < MIN_QUAL_BAD_OVERLAP && secondReadQuals[iii-firstReadStop] < MIN_QUAL_BAD_OVERLAP ) {
                return overlappingPair; // both reads have low qual bases in the overlap region so don't merge them because don't know what is going on
            }
            bases[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadBases[iii] : secondReadBases[iii-firstReadStop] );
            quals[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadQuals[iii] : secondReadQuals[iii-firstReadStop] );
        }
        for(int iii = firstRead.getReadLength(); iii < numBases; iii++) {
            bases[iii] = secondReadBases[iii-firstReadStop];
            quals[iii] = secondReadQuals[iii-firstReadStop];
        }

        final GATKSAMRecord returnRead = new GATKSAMRecord( firstRead.getHeader() );
        returnRead.setAlignmentStart( firstRead.getUnclippedStart() );
        returnRead.setReadBases( bases );
        returnRead.setBaseQualities( quals );
        returnRead.setReadGroup( firstRead.getReadGroup() );
        returnRead.setReferenceName( firstRead.getReferenceName() );
        returnRead.setReadName( firstRead.getReadName() );
        final CigarElement c = new CigarElement(bases.length, CigarOperator.M);
        final ArrayList<CigarElement> cList = new ArrayList<CigarElement>();
        cList.add(c);
        returnRead.setCigar( new Cigar( cList ));
        returnRead.setMappingQuality( firstRead.getMappingQuality() );

        if( firstRead.hasBaseIndelQualities() || secondRead.hasBaseIndelQualities() ) {
            final byte[] firstReadInsertionQuals = firstRead.getBaseInsertionQualities();
            final byte[] firstReadDeletionQuals = firstRead.getBaseDeletionQualities();
            final byte[] secondReadInsertionQuals = secondRead.getBaseInsertionQualities();
            final byte[] secondReadDeletionQuals = secondRead.getBaseDeletionQualities();
            for(int iii = 0; iii < firstReadStop; iii++) {
                insertionQuals[iii] = firstReadInsertionQuals[iii];
                deletionQuals[iii] = firstReadDeletionQuals[iii];
            }
            for(int iii = firstReadStop; iii < firstRead.getReadLength(); iii++) {
                insertionQuals[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadInsertionQuals[iii] : secondReadInsertionQuals[iii-firstReadStop] ); // Purposefully checking the highest *base* quality score
                deletionQuals[iii] = ( firstReadQuals[iii] > secondReadQuals[iii-firstReadStop] ? firstReadDeletionQuals[iii] : secondReadDeletionQuals[iii-firstReadStop] ); // Purposefully checking the highest *base* quality score
            }
            for(int iii = firstRead.getReadLength(); iii < numBases; iii++) {
                insertionQuals[iii] = secondReadInsertionQuals[iii-firstReadStop];
                deletionQuals[iii] = secondReadDeletionQuals[iii-firstReadStop];
            }
            returnRead.setBaseQualities( insertionQuals, EventType.BASE_INSERTION );
            returnRead.setBaseQualities( deletionQuals, EventType.BASE_DELETION );
        }

        final ArrayList<GATKSAMRecord> returnList = new ArrayList<GATKSAMRecord>();
        returnList.add(returnRead);
        return returnList;
    }
}
