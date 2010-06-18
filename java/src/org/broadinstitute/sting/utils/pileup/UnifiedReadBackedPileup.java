package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.gatk.iterators.IterableIterator;
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
public class UnifiedReadBackedPileup implements ReadBackedPileup {
    private GenomeLoc loc = null;
    private List<PileupElement> pileup = null;
    
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
    public UnifiedReadBackedPileup(GenomeLoc loc, List<SAMRecord> reads, List<Integer> offsets ) {
        this(loc, readsOffsets2Pileup(reads, offsets));
    }

    public UnifiedReadBackedPileup(GenomeLoc loc, List<SAMRecord> reads, int offset ) {
        this(loc, readsOffsets2Pileup(reads, offset));
    }

    /**
     * Create a new version of a read backed pileup at loc without any aligned reads
     *
     */
     public UnifiedReadBackedPileup(GenomeLoc loc ) {
        this(loc, new ArrayList<PileupElement>(0));
    }

    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.  This lower level constructure assumes pileup is well-formed and merely keeps a
     * pointer to pileup.  Don't go changing the data in pileup.
     *
     */
     public UnifiedReadBackedPileup(GenomeLoc loc, List<PileupElement> pileup ) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in ReadBackedPileup");
        if ( pileup == null ) throw new StingException("Illegal null pileup in ReadBackedPileup");

        this.loc = loc;
        this.pileup = pileup;

        calculatedCachedData();
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc
     * @param pileup
     */
    public UnifiedReadBackedPileup(GenomeLoc loc, List<PileupElement> pileup, int size, int nDeletions, int nMQ0Reads ) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in UnifiedReadBackedPileup");
        if ( pileup == null ) throw new StingException("Illegal null pileup in UnifiedReadBackedPileup");

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
        if ( reads == null ) throw new StingException("Illegal null read list in UnifiedReadBackedPileup");
        if ( offsets == null ) throw new StingException("Illegal null offsets list in UnifiedReadBackedPileup");
        if ( reads.size() != offsets.size() ) throw new StingException("Reads and offset lists have different sizes!");

        ArrayList<PileupElement> pileup = new ArrayList<PileupElement>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            pileup.add( new PileupElement( reads.get(i), offsets.get(i) ) );
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
    private static ArrayList<PileupElement> readsOffsets2Pileup(List<SAMRecord> reads, int offset ) {
        if ( reads == null ) throw new StingException("Illegal null read list in UnifiedReadBackedPileup");
        if ( offset < 0 ) throw new StingException("Illegal offset < 0 UnifiedReadBackedPileup");

        ArrayList<PileupElement> pileup = new ArrayList<PileupElement>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            pileup.add( new PileupElement( reads.get(i), offset ) );
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
    @Override
    public ReadBackedPileup getPileupWithoutDeletions() {
        if ( getNumberOfDeletions() > 0 ) {
            ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();

            for ( PileupElement p : pileup ) {
                if ( !p.isDeletion() ) {
                    filteredPileup.add(p);
                }
            }
            return new UnifiedReadBackedPileup(loc, filteredPileup);
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
        Map<String, PileupElement> filteredPileup = new HashMap<String, PileupElement>();

        for ( PileupElement p : pileup ) {
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

        return new UnifiedReadBackedPileup(loc, new ArrayList<PileupElement>(filteredPileup.values()));
    }

    /**
     * Returns a new ReadBackedPileup that is free of mapping quality zero reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no MQ0 reads in the pileup.
     *
     * @return
     */
    @Override
    public ReadBackedPileup getPileupWithoutMappingQualityZeroReads() {

        if ( getNumberOfMappingQualityZeroReads() > 0 ) {
            ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();
            for ( PileupElement p : pileup ) {
                if ( p.getRead().getMappingQuality() > 0 ) {
                    filteredPileup.add(p);
                }
            }
            return new UnifiedReadBackedPileup(loc, filteredPileup);
        } else {
            return this;
        }
    }

    /** Returns subset of this pileup that contains only bases with quality >= minBaseQ, coming from
     * reads with mapping qualities >= minMapQ. This method allocates and returns a new instance of ReadBackedPileup.
     * @param minBaseQ
     * @param minMapQ
     * @return
     */
    @Override
    public ReadBackedPileup getBaseAndMappingFilteredPileup( int minBaseQ, int minMapQ ) {
        ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();

        for ( PileupElement p : pileup ) {
            if ( p.getRead().getMappingQuality() >= minMapQ && (p.isDeletion() || p.getQual() >= minBaseQ) ) {
                filteredPileup.add(p);
            }
        }

        return new UnifiedReadBackedPileup(loc, filteredPileup);
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
            if ( positions.add(generator.nextInt(pileup.size())) )
                i++;
        }

        Iterator positionIter = positions.iterator();
        ArrayList<PileupElement> downsampledPileup = new ArrayList<PileupElement>();

        while ( positionIter.hasNext() ) {
            int nextReadToKeep = (Integer)positionIter.next();
            downsampledPileup.add(pileup.get(nextReadToKeep));
        }

        return new UnifiedReadBackedPileup(getLocation(), downsampledPileup);
    }

    @Override
    public Collection<String> getSamples() {
        Collection<String> sampleNames = new HashSet<String>();
        for(PileupElement p: this) {
            SAMRecord read = p.getRead();
            String sampleName = read.getReadGroup() != null ? read.getReadGroup().getSample() : null;
            sampleNames.add(sampleName);
        }
        return sampleNames;
    }

    @Override
    public ReadBackedPileup getPileupForSample(String sampleName) {
        List<PileupElement> filteredPileup = new ArrayList<PileupElement>();
        for(PileupElement p: this) {
            SAMRecord read = p.getRead();
            if(sampleName != null) {
                if(read.getReadGroup() != null && sampleName.equals(read.getReadGroup().getSample()))
                    filteredPileup.add(p);
            }
            else {
                if(read.getReadGroup() == null || read.getReadGroup().getSample() == null)
                    filteredPileup.add(p);
            }
        }
        return filteredPileup.size()>0 ? new UnifiedReadBackedPileup(loc,filteredPileup) : null;
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
    // todo -- why is this public?
    @Override
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
    @Override
    public boolean hasSecondaryBases() {
        for ( PileupElement pile : this ) {
             // skip deletion sites
             if ( ! pile.isDeletion() && BaseUtils.isRegularBase((char)pile.getSecondBase()) )
                 return true;
        }

        return false;
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


    //
    // Private functions for printing pileups
    //
    private String getMappingQualsString() {
        return quals2String(getMappingQuals());
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
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(getLocation());
        s.append(": ");

        for ( PileupElement p : this ) {
            s.append((char)p.getBase());
        }

        return s.toString();
    }
}
