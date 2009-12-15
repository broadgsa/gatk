package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.iterators.IterableIterator;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ExtendedPileupElement;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.*;

import net.sf.samtools.SAMRecord;

/**
 * Version two file implementing pileups of bases in reads at a locus.
 *
 * @author Mark DePristo
 */
public class ReadBackedPileup implements Iterable<PileupElement> {
    private GenomeLoc loc = null;
    private ArrayList<PileupElement> pileup = null;
    
    private int size = 0;                   // cached value of the size of the pileup
    private int nDeletions = 0;             // cached value of the number of deletions
    private int nMQ0Reads = 0;              // cached value of the number of MQ0 reads

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
    public ReadBackedPileup(GenomeLoc loc, List<SAMRecord> reads, List<Integer> offsets ) {
        this(loc, readsOffsets2Pileup(reads, offsets));
    }


    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.  This lower level constructure assumes pileup is well-formed and merely keeps a
     * pointer to pileup.  Don't go changing the data in pileup.
     *
     */
     public ReadBackedPileup(GenomeLoc loc, ArrayList<PileupElement> pileup ) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in ReadBackedPileup2");
        if ( pileup == null ) throw new StingException("Illegal null pileup in ReadBackedPileup2");

        this.loc = loc;
        this.pileup = pileup;

        calculatedCachedData();
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc
     * @param pileup
     */
    public ReadBackedPileup(GenomeLoc loc, ArrayList<PileupElement> pileup, int size, int nDeletions, int nMQ0Reads ) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in ReadBackedPileup2");
        if ( pileup == null ) throw new StingException("Illegal null pileup in ReadBackedPileup2");

        this.loc = loc;
        this.pileup = pileup;
        this.size = size;
        this.nDeletions = nDeletions;
        this.nMQ0Reads = nMQ0Reads;
   }


    /**
     * Calculate cached sizes, nDeletion, and base counts for the pileup.  This calculation is done upfront,
     * so you pay the cost at the start, but it's more efficient to do this rather than pay the cost of calling
     * sizes, nDeletion, etc. over and over potentially.
     */
    private void calculatedCachedData() {
        size = 0;
        nDeletions = 0;
        nMQ0Reads = 0;

        for ( PileupElement p : this ) {
            size++;
            if ( p.isDeletion() ) {
                nDeletions++;
            }
            if ( p.getRead().getMappingQuality() == 0 ) {
                nMQ0Reads++;
            }
        }
    }


    /**
     * Helper routine for converting reads and offset lists to a PileupElement list.
     *
     * @param reads
     * @param offsets
     * @return
     */
    private static ArrayList<PileupElement> readsOffsets2Pileup(List<SAMRecord> reads, List<Integer> offsets ) {
        if ( reads == null ) throw new StingException("Illegal null read list in ReadBackedPileup2");
        if ( offsets == null ) throw new StingException("Illegal null offsets list in ReadBackedPileup2");
        if ( reads.size() != offsets.size() ) throw new StingException("Reads and offset lists have different sizes!");

        ArrayList<PileupElement> pileup = new ArrayList<PileupElement>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            pileup.add( new PileupElement( reads.get(i), offsets.get(i) ) );
        }

        return pileup;
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
    public ReadBackedPileup getPileupWithoutDeletions() {
        if ( getNumberOfDeletions() > 0 ) {
            ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();

            for ( PileupElement p : pileup ) {
                if ( !p.isDeletion() ) {
                    filteredPileup.add(p);
                }
            }
            return new ReadBackedPileup(loc, filteredPileup);
        } else {
            return this;
        }
    }

    /**
     * Returns a new ReadBackedPileup that is free of mapping quality zero reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no MQ0 reads in the pileup.
     *
     * @return
     */
    public ReadBackedPileup getPileupWithoutMappingQualityZeroReads() {

        if ( getNumberOfMappingQualityZeroReads() > 0 ) {
            ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();
            for ( PileupElement p : pileup ) {
                if ( p.getRead().getMappingQuality() > 0 ) {
                    filteredPileup.add(p);
                }
            }
            return new ReadBackedPileup(loc, filteredPileup);
        } else {
            return this;
        }
    }

    public ReadBackedPileup getBaseAndMappingFilteredPileup( int minBaseQ, int minMapQ ) {
        ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();

        for ( PileupElement p : pileup ) {
            if ( p.getRead().getMappingQuality() >= minMapQ && (p.isDeletion() || p.getQual() >= minBaseQ) ) {
                filteredPileup.add(p);
            }
        }

        return new ReadBackedPileup(loc, filteredPileup);
    }
    /**
     * Returns a pileup randomly downsampled to the desiredCoverage.
     *
     * @param desiredCoverage
     * @return
     */
    public ReadBackedPileup getDownsampledPileup(int desiredCoverage) {
        if ( size() <= desiredCoverage )
            return this;

        // randomly choose numbers corresponding to positions in the reads list
        Random generator = new Random();
        TreeSet<Integer> positions = new TreeSet<Integer>();
        for ( int i = 0; i < desiredCoverage; /* no update */ ) {
            if ( positions.add(generator.nextInt(pileup.size())) )
                i++;
        }

        Iterator positionIter = positions.iterator();
        ArrayList<PileupElement> downsampledPileup = new ArrayList<PileupElement>();

        while ( positionIter.hasNext() ) {
            int nextReadToKeep = (Integer)positionIter.next();
            downsampledPileup.add(pileup.get(nextReadToKeep));
        }

        return new ReadBackedPileup(getLocation(), downsampledPileup);
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
    public Iterator<PileupElement> iterator() {
        return pileup.iterator();
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

    // todo -- reimplement efficiently
    public IterableIterator<ExtendedPileupElement> extendedForeachIterator() {
        ArrayList<ExtendedPileupElement> x = new ArrayList<ExtendedPileupElement>(size());
        int i = 0;
        for ( PileupElement pile : this ) {
            x.add(new ExtendedPileupElement(pile.getRead(), pile.getOffset(), i++, this)); 
        }

        return new IterableIterator<ExtendedPileupElement>(x.iterator());
    }

    /**
     * Simple useful routine to count the number of deletion bases in this pileup
     *
     * @return
     */
    public int getNumberOfDeletions() {
        return nDeletions;
    }

    public int getNumberOfMappingQualityZeroReads() {
        return nMQ0Reads;
    }

//    public int getNumberOfDeletions() {
//        int n = 0;
//
//        for ( int i = 0; i < size(); i++ ) {
//            if ( getOffsets().get(i) != -1 ) { n++; }
//        }
//
//        return n;
//    }

    /**
     * @return the number of elements in this pileup
     */
    public int size() {
        return size;
    }

    /**
     * @return the location of this pileup
     */
    public GenomeLoc getLocation() {
        return loc;
    }

    /**
     * Get counts of A, C, G, T in order, which returns a int[4] vector with counts according
     * to BaseUtils.simpleBaseToBaseIndex for each base.
     *
     * @return
     */
    public int[] getBaseCounts() {
        int[] counts = new int[4];
        for ( PileupElement pile : this ) {
            // skip deletion sites
            if ( ! pile.isDeletion() ) {
                int index = BaseUtils.simpleBaseToBaseIndex((char)pile.getBase());
                if (index != -1)
                    counts[index]++;
            }
        }

        return counts;
    }

    /**
     * Somewhat expensive routine that returns true if any base in the pileup has secondary bases annotated
     * @return
     */
    public boolean hasSecondaryBases() {
        for ( PileupElement pile : this ) {
             // skip deletion sites
             if ( ! pile.isDeletion() && BaseUtils.isRegularBase((char)pile.getSecondBase()) )
                 return true;
        }

        return false;
    }

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
    public List<SAMRecord> getReads() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>(size());
        for ( PileupElement pile : this ) { reads.add(pile.getRead()); }
        return reads;
    }

    /**
     * Returns a list of the offsets in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    public List<Integer> getOffsets() {
        List<Integer> offsets = new ArrayList<Integer>(size());
        for ( PileupElement pile : this ) { offsets.add(pile.getOffset()); }
        return offsets;
    }

    /**
     * Returns an array of the bases in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    public byte[] getBases() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getBase(); }
        return v;
    }

    /**
     * Returns an array of the secondary bases in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    public byte[] getSecondaryBases() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getSecondBase(); }
        return v;
    }

     /**
     * Returns an array of the quals in this pileup. Note this call costs O(n) and allocates fresh array each time
     * @return
     */
     public byte[] getQuals() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getQual(); }
        return v;
    }

    /**
     * Get an array of the mapping qualities
     * @return
     */
    public byte[] getMappingQuals() {
       byte[] v = new byte[size()];
       for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = (byte)pile.getRead().getMappingQuality(); }
       return v;
   }


    //
    // Private functions for printing pileups
    //
    private String getMappingQualsString() {
        return quals2String(getMappingQuals());
    }

    private static String quals2String( byte[] quals ) {
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
}
