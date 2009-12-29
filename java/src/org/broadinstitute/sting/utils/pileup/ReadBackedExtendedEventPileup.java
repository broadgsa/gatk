package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;

import java.util.*;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Dec 29, 2009
 * Time: 12:25:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadBackedExtendedEventPileup implements Iterable<ExtendedEventPileupElement>  {
    private GenomeLoc loc = null;
    private ArrayList<ExtendedEventPileupElement> pileup = null;

    private int size = 0;                   // cached value of the size of the pileup
    private int nDeletions = 0;             // cached value of the number of deletions
    private int nInsertions = 0;
    private int nMQ0Reads = 0;              // cached value of the number of MQ0 reads

    /**
     * Create a new version of a read backed pileup at loc without any aligned reads
     *
     */
     public ReadBackedExtendedEventPileup(GenomeLoc loc ) {
        this(loc, new ArrayList<ExtendedEventPileupElement>(0));
    }

    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.  This lower level constructure assumes pileup is well-formed and merely keeps a
     * pointer to pileup.  Don't go changing the data in pileup.
     *
     */
     public ReadBackedExtendedEventPileup(GenomeLoc loc, ArrayList<ExtendedEventPileupElement> pileup ) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in ReadBackedExtendedEventPileup");
        if ( pileup == null ) throw new StingException("Illegal null pileup in ReadBackedExtendedEventPileup");

        this.loc = loc;
        this.pileup = pileup;

        calculatedCachedData();
    }

    /**
     * Optimization of above constructor where all of the cached data is provided
     * @param loc
     * @param pileup
     */
    public ReadBackedExtendedEventPileup(GenomeLoc loc, ArrayList<ExtendedEventPileupElement> pileup, int size, int nInsertions, int nDeletions, int nMQ0Reads ) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in ReadBackedExtendedEventPileup");
        if ( pileup == null ) throw new StingException("Illegal null pileup in ReadBackedExtendedEventPileup");

        this.loc = loc;
        this.pileup = pileup;
        this.size = size;
        this.nDeletions = nDeletions;
        this.nInsertions = nInsertions;
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
        nInsertions = 0;
        nMQ0Reads = 0;

        for ( ExtendedEventPileupElement p : this ) {
            size++;
            if ( p.isDeletion() ) {
                nDeletions++;
            } else {
                if ( p.isInsertion() ) nInsertions++;
            }

            if ( p.getRead().getMappingQuality() == 0 ) {
                nMQ0Reads++;
            }
        }
    }


    // --------------------------------------------------------
    //
    // Special 'constructors'
    //
    // --------------------------------------------------------


    /**
     * Returns a new ReadBackedPileup that is free of mapping quality zero reads in this pileup.  Note that this
     * does not copy the data, so both ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
     * of the pileup (just returns this) if there are no MQ0 reads in the pileup.
     *
     * @return
     */
    public ReadBackedExtendedEventPileup getPileupWithoutMappingQualityZeroReads() {

        if ( getNumberOfMappingQualityZeroReads() > 0 ) {
            ArrayList<ExtendedEventPileupElement> filteredPileup = new ArrayList<ExtendedEventPileupElement>();
            for ( ExtendedEventPileupElement p : pileup ) {
                if ( p.getRead().getMappingQuality() > 0 ) {
                    filteredPileup.add(p);
                }
            }
            return new ReadBackedExtendedEventPileup(loc, filteredPileup);
        } else {
            return this;
        }
    }

    public ReadBackedExtendedEventPileup getMappingFilteredPileup( int minMapQ ) {
        ArrayList<ExtendedEventPileupElement> filteredPileup = new ArrayList<ExtendedEventPileupElement>();

        for ( ExtendedEventPileupElement p : pileup ) {
            if ( p.getRead().getMappingQuality() >= minMapQ  ) {
                filteredPileup.add(p);
            }
        }

        return new ReadBackedExtendedEventPileup(loc, filteredPileup);
    }

    /**
     * Returns a pileup randomly downsampled to the desiredCoverage.
     *
     * @param desiredCoverage
     * @return
     */
    public ReadBackedExtendedEventPileup getDownsampledPileup(int desiredCoverage) {
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
        ArrayList<ExtendedEventPileupElement> downsampledPileup = new ArrayList<ExtendedEventPileupElement>();

        while ( positionIter.hasNext() ) {
            int nextReadToKeep = (Integer)positionIter.next();
            downsampledPileup.add(pileup.get(nextReadToKeep));
        }

        return new ReadBackedExtendedEventPileup(getLocation(), downsampledPileup);
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
    public Iterator<ExtendedEventPileupElement> iterator() {
        return pileup.iterator();
    }


    /**
     * Returns the number of deletion events in this pileup
     *
     * @return
     */
    public int getNumberOfDeletions() {
        return nDeletions;
    }

    /**
     * Returns the number of insertion events in this pileup
     *
     * @return
     */
    public int getNumberOfInsertions() {
        return nInsertions;
    }

    /**
     * Returns the number of mapping quality zero reads in this pileup.
     * @return
     */
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



    public String getShortPileupString() {
        // In the pileup format, each extended event line has genomic position (chromosome name and offset),
        // reference "base" (always set to "E" for E(xtended)), then 'I','D', or '.' symbol for each read representing
        // insertion, deletion or no-event, respectively.
        return String.format("%s %s E %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                new String(getEvents()) );
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
        for ( ExtendedEventPileupElement pile : this ) { reads.add(pile.getRead()); }
        return reads;
    }

    /**
     * Returns a list of the offsets in this pileup. Note this call costs O(n) and allocates fresh lists each time
     * @return
     */
    public List<Integer> getOffsets() {
        List<Integer> offsets = new ArrayList<Integer>(size());
        for ( ExtendedEventPileupElement pile : this ) { offsets.add(pile.getOffset()); }
        return offsets;
    }

    /**
     * Returns an array of the events in this pileup ('I', 'D', or '.'). Note this call costs O(n) and allocates fresh array each time
     * @return
     */
    public byte[] getEvents() {
        byte[] v = new byte[size()];
        int i = 0;
        for ( ExtendedEventPileupElement e : this ) {
            switch ( e.getType() ) {
                case INSERTION: v[i] = 'I'; break;
                case DELETION: v[i] = 'D'; break;
                case NOEVENT: v[i] = '.'; break;
                default: throw new StingException("Unknown event type encountered: "+e.getType());
            }
            i++;
        }
        return v;
    }

    public List<Pair<String,Integer>> getEventStringsWithCounts() {
        Map<String, Integer> events = new HashMap<String,Integer>();

        for ( ExtendedEventPileupElement e : this ) {
            Integer cnt;
            String indel = null;
            switch ( e.getType() ) {
                case INSERTION:
                    indel = "+"+e.getEventBases();
                    break;
                case DELETION:
                    indel = Integer.toString(e.getEventLength())+"D";
                    break;
                case NOEVENT: continue;
                default: throw new StingException("Unknown event type encountered: "+e.getType());
            }

            cnt = events.get(indel);
            if ( cnt == null ) events.put(indel,1);
            else events.put(indel,cnt.intValue()+1);
        }

        List<Pair<String,Integer>> eventList = new ArrayList<Pair<String,Integer>>(events.size());
        for ( Map.Entry<String,Integer> m : events.entrySet() ) {
            eventList.add( new Pair<String,Integer>(m.getKey(),m.getValue()));
        }
        return eventList;
    }


    /**
     * Get an array of the mapping qualities
     * @return
     */
    public byte[] getMappingQuals() {
       byte[] v = new byte[size()];
       int i = 0;
       for ( ExtendedEventPileupElement e : this ) { v[i++] = (byte)e.getRead().getMappingQuality(); }
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

}
