package org.broadinstitute.sting.gatk.iterators;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.util.RuntimeIOException;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.io.IOException;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.FastaSequenceFile2;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:45:01 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceIterator implements Iterator<ReferenceIterator> {

    // enable debugging output?
    private final boolean DEBUG = false;

    // The reference sequence file generator
    private FastaSequenceFile2 refFile;

    private ReferenceSequence currentContig = null;
    //private ReferenceSequence nextContig = null;    
    private long offset = -1;

    public ReferenceIterator( FastaSequenceFile2 refFile ) {
        this.refFile = refFile;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Accessing data
    //
    // --------------------------------------------------------------------------------------------------------------
    public byte getBaseAsByte() { return currentContig.getBases()[(int)offset]; }
    public String getBaseAsString() { 
	assert offset > -1 : currentContig.getName() + " index is " + offset;
	//assert offset < currentContig.getBases().();

	return StringUtil.bytesToString(currentContig.getBases(), (int)offset, 1);
    }
    public char getBaseAsChar() { return getBaseAsString().charAt(0); }
    public ReferenceSequence getCurrentContig() { return currentContig; }
    public long getPosition() { return offset + 1; }
    public GenomeLoc getLocation() { return new GenomeLoc( getCurrentContig().getName(), getPosition() ); }
    
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
            return readNextContig();
        }
    }

    public ReferenceIterator next() {
        if ( currentContig != null ) {
            if ( DEBUG ) System.out.printf("  -> %s:%d %d%n", currentContig.getName(), offset, currentContig.length());
        }
        offset++;  // move on to the next position

        if ( currentContig == null || offset >= currentContig.length() ) {
            // We need to update the contig
            if ( readNextContig() ){
                // We sucessfully loaded the next contig, recursively call next
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

    /**
     * Helper routine that doesn't move the contigs around, it just checks that everything is kosher in the seek
     * within this chromosome
     * @param seekContigName name for printing pursues, asserted to be the current contig name
     * @param seekOffset where we want to be in this contig
     * @return this setup to be at seekoffset within seekContigName 
     */
    private ReferenceIterator seekForwardOffsetOnSameContig(final String seekContigName, final long seekOffset) {
        assert seekContigName.equals(currentContig.getName()) : String.format("only works on this contig, but the current %s and sought %s contigs are different!", currentContig.getName(), seekContigName);

        // we're somewhere on this contig
        if ( seekOffset < offset || seekOffset >= currentContig.length() ) {
            // bad boy -- can't go backward safely or just beyond the contig length
            throw new IllegalArgumentException(String.format("Invalid seek to %s from %s, which is usually due to out of order reads%n",
                    new GenomeLoc(currentContig.getName(), seekOffset), new GenomeLoc(currentContig.getName(), offset)));
        }
        else {
            offset = seekOffset - 1;
            return next();
        }
    }


    private ReferenceIterator seekForwardOffset(final String seekContigName, final long seekOffset) {
        assert seekContigName != null : "seekContigName is null";
        assert seekOffset >= 0 : "seekOffset < 0: " + seekOffset;

        // jumps us forward in the sequence to the contig / pos
        if ( currentContig == null )
            next();

        if ( DEBUG ) System.out.printf("  -> Seeking to %s %d from %s %d%n", seekContigName, seekOffset, currentContig.getName(), offset);

        int cmpContigs = GenomeLoc.compareContigs(seekContigName, currentContig.getName());

        if ( cmpContigs == -1 && false ) {  // todo: fixed
            // The contig we are looking for is before the currentContig -- it's an error
            throw new IllegalArgumentException(String.format("Invalid seek to %s from %s, which is usually due to out of order reads%n",
                    new GenomeLoc(currentContig.getName(), seekOffset), new GenomeLoc(currentContig.getName(), offset)));
        }
        else if ( cmpContigs == 1 ) {
            // we need to jump forward
            if ( DEBUG ) System.out.printf("  -> Seeking in the fasta file to %s from %s%n", seekContigName, currentContig.getName());

            if ( ! refFile.seekToContig(seekContigName) ) { // ok, do the seek
                // a false result indicates a failure, throw a somewhat cryptic call
                throw new RuntimeIOException(String.format("Unexpected seek failure from %s from %s%n",
                        new GenomeLoc(currentContig.getName(), seekOffset), new GenomeLoc(currentContig.getName(), offset)));
            }

            readNextContig(); // since we haven't failed, we just read in the next contig (which is seekContigName)
        }

        // at this point, the current contig is seekContigName, so just do a bit more error checking and be done
        return seekForwardOffsetOnSameContig( seekContigName, seekOffset );
    }


    // --------------------------------------------------------------------------------------------------------------
    //
    // Interal state manipulation
    //
    // --------------------------------------------------------------------------------------------------------------
    protected boolean readNextContig() {
        // returns true if we had another contig to load
        currentContig = refFile.nextSequence();
        offset = -1;
        return currentContig != null;
    }

    /**
     * Simple forwarding method to the refFile itself
     * 
     * @return
     */
    public String nextContigName() {
        return this.refFile.getNextContigName();
    }

}
