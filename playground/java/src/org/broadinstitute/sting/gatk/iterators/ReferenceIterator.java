package org.broadinstitute.sting.gatk.iterators;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

import java.util.Iterator;
import java.util.NoSuchElementException;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:45:01 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceIterator implements Iterator<ReferenceIterator> {

    // The reference sequence file generator
    private ReferenceSequenceFile refFile;

    private ReferenceSequence currentContig = null;
    private ReferenceSequence nextContig = null;    
    private long offset = -1;

    public ReferenceIterator( ReferenceSequenceFile refFile ) {
        this.refFile = refFile;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Accessing data
    //
    // --------------------------------------------------------------------------------------------------------------
    public byte getBaseAsByte() { return currentContig.getBases()[(int)offset]; }
    public String getBaseAsString() { return StringUtil.bytesToString(currentContig.getBases(), (int)offset, 1); }
    public char getBaseAsChar() { return getBaseAsString().charAt(0); }
    public ReferenceSequence getCurrentContig() { return currentContig; }
    public long getPosition() { return offset + 1; }
    
    // --------------------------------------------------------------------------------------------------------------
    //
    // Iterator routines
    //
    // --------------------------------------------------------------------------------------------------------------
    public boolean hasNext() {
        if ( currentContig == null || offset + 1 < currentContig.length() ) {
            return true;
        }
        else {
            return loadNextContig();
        }
    }

    public ReferenceIterator next() {
        if ( currentContig != null ) {
            //System.out.printf("  -> %s:%d %d%n", currentContig.getName(), offset, currentContig.length());
        }
        offset++;  // move on to the next position

        if ( currentContig == null || offset >= currentContig.length() ) {
            // We need to update the contig
            //System.out.printf("  -> Updating length%n");
            if ( nextContig != null ) {
                // We've already loaded the next contig, swap it in, and recursively call next
                swapNextContig();
                return next();
            }
            else if ( loadNextContig() ){
                // We sucessfully loaded the next contig, recursively call next
                offset = -1;
                return next();
            }
            else {
                throw new NoSuchElementException();
            }
        }
        else {
            // We're good to go -- we're in the current contig
            return this;
        }
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    
    // --------------------------------------------------------------------------------------------------------------
    //
    // Jumping forward
    //
    // --------------------------------------------------------------------------------------------------------------
    public ReferenceIterator seekForward(final GenomeLoc loc) {
        assert loc != null : "seekForward location is null";
        
        return seekForwardOffset(loc.getContig(), loc.getStart() - 1);
    }

    public ReferenceIterator seekForward(final String contigName, final long pos) {
        return seekForwardOffset(contigName, pos - 1);
    }

    private ReferenceIterator seekForwardOffset(final String contigName, final long seekOffset) {
        assert contigName != null : "contigName is null";
        assert seekOffset >= 0 : "seekOffset < 0: " + seekOffset;

        // jumps us forward in the sequence to the contig / pos
        if ( currentContig == null )
            next();

        //System.out.printf("  -> Seeking to %s %d from %s %d%n", contigName, seekOffset, currentContig.getName(), offset);
        if ( contigName.equals(currentContig.getName()) ) {
            // we're somewhere on this contig
            if ( seekOffset < offset || seekOffset >= currentContig.length() ) {
                // bad boy -- can't go backward safely or just beyond the contig length
                throw new IllegalArgumentException(String.format("Invalid seek to %s from %s, which is usually due to out of order reads%n",
                        new GenomeLoc(currentContig.getName(), seekOffset), new GenomeLoc(currentContig.getName(), offset)));
                //return null;
            }
            else {
                offset = seekOffset - 1;
                return next();
            }
        }
        else {
            while (true) {
                //System.out.printf("Seeking to contig %s, cur=%s, next=%s%n", contigName, currentContig.getName(),
                //        nextContig != null ? nextContig.getName() : "not loaded yet");
                // go searching through the reference
                if ( ! loadNextContig() ) {
                    // never found anything
                    return null;
                }
                else if ( nextContig.getName().equals(contigName) ) {
                    swapNextContig();
                    return seekForward(contigName, seekOffset);
                }
            }
        }
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // Interal state manipulation
    //
    // --------------------------------------------------------------------------------------------------------------
    protected boolean loadNextContig() {
        // returns true if we had another contig to load
        nextContig = refFile.nextSequence();
        return nextContig != null;
    }

    protected void swapNextContig() {
        currentContig = nextContig;
        nextContig = null;
        offset = -1;
    }
}
