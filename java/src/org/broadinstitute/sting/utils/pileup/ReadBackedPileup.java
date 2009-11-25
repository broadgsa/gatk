package org.broadinstitute.sting.utils.pileup;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.iterators.IterableIterator;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ExtendedPileupElement;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

import net.sf.samtools.SAMRecord;

/**
 * Version two file implementing pileups of bases in reads at a locus.
 *
 * @author Mark DePristo
 */
public class ReadBackedPileup implements Iterable<PileupElement> {
    private GenomeLoc loc = null;
    private char ref = 0;
    private ArrayList<PileupElement> pileup = null;

    public ReadBackedPileup(char ref, AlignmentContext context ) {
        this(context.getLocation(), ref, context.getReads(), context.getOffsets());
    }

    public ReadBackedPileup(GenomeLoc loc, char ref, List<SAMRecord> reads, List<Integer> offsets ) {
        this(loc, ref, readsOffsets2Pileup(reads, offsets));
    }

    public ReadBackedPileup(GenomeLoc loc, char ref, ArrayList<PileupElement> pileup ) {
        if ( loc == null ) throw new StingException("Illegal null genomeloc in ReadBackedPileup2");
        if ( pileup == null ) throw new StingException("Illegal null pileup in ReadBackedPileup2");

        this.loc = loc;
        this.ref = ref;
        this.pileup = pileup;
    }

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

    //
    // iterators
    //
    public Iterator<PileupElement> iterator() {
        return pileup.iterator();
    }

    // todo -- reimplement efficiently
    public IterableIterator<ExtendedPileupElement> extendedForeachIterator() {
        ArrayList<ExtendedPileupElement> x = new ArrayList<ExtendedPileupElement>(size());
        int i = 0;
        for ( PileupElement pile : this ) {
            x.add(new ExtendedPileupElement(pile.getRead(), pile.getOffset(), i++, this)); 
        }

        return new IterableIterator<ExtendedPileupElement>(x.iterator());
    }


    public ReadBackedPileup getPileupWithoutDeletions() {
        // todo -- fixme
        if ( getNumberOfDeletions() > 0 ) { // todo -- remember number of deletions
            List<SAMRecord> newReads = new ArrayList<SAMRecord>();
            List<Integer> newOffsets = new ArrayList<Integer>();

            for ( int i = 0; i < size(); i++ ) {
                if ( getOffsets().get(i) != -1 ) {
                    newReads.add(getReads().get(i));
                    newOffsets.add(getOffsets().get(i));
                }
            }
            return new ReadBackedPileup(loc, ref, newReads, newOffsets);
        } else {
            return this;
        }
    }

    public int getNumberOfDeletions() {
        int n = 0;

        for ( int i = 0; i < size(); i++ ) {
            if ( getOffsets().get(i) != -1 ) { n++; }
        }

        return n;
    }

    // todo -- optimize me
    public int size() {
        return pileup.size();
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public char getRef() {
        return ref;
    }

    public List<SAMRecord> getReads() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>(size());
        for ( PileupElement pile : this ) { reads.add(pile.getRead()); }
        return reads;
    }

    public List<Integer> getOffsets() {
        List<Integer> offsets = new ArrayList<Integer>(size());
        for ( PileupElement pile : this ) { offsets.add(pile.getOffset()); }
        return offsets;
    }

    public byte[] getBases() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getBase(); }
        return v;
    }

    public byte[] getSecondaryBases() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getSecondBase(); }
        return v;
    }

    public byte[] getQuals() {
        byte[] v = new byte[size()];
        for ( ExtendedPileupElement pile : this.extendedForeachIterator() ) { v[pile.getPileupOffset()] = pile.getQual(); }
        return v;
    }

    public int[] getBaseCounts() {
        int[] counts = new int[4];
        for ( PileupElement pile : this ) {
            // skip deletion sites
            if ( ! pile.isDeletion() ) {
                char base = Character.toUpperCase((char)(pile.getBase()));
                if (BaseUtils.simpleBaseToBaseIndex(base) == -1)
                    continue;
                counts[BaseUtils.simpleBaseToBaseIndex(base)]++;
            }
        }

        return counts;
    }

    public boolean hasSecondaryBases() {
        for ( PileupElement pile : this ) {
             // skip deletion sites
             if ( ! pile.isDeletion() && BaseUtils.isRegularBase((char)pile.getSecondBase()) )
                 return true;
        }

        return false;
    }

    public String getPileupString(boolean qualsAsInts) {
        // In the pileup format, each line represents a genomic position, consisting of chromosome name,
        // coordinate, reference base, read bases, read qualities and alignment mapping qualities.

        //return String.format("%s %s %s %s", getLocation(), getRef(), getBases(), getQuals());
        return String.format("%s %s %s %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                getRef(),                                               // reference base
                new String(getBases()));
    }
}
