package org.broadinstitute.sting.utils.fasta;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.io.IoUtil;

import java.io.*;

import net.sf.samtools.SAMTextHeaderCodec;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.AsciiLineReader;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.util.RuntimeIOException;
import org.apache.log4j.Logger;

/**
 * Implementation of ReferenceSequenceFile for reading from FASTA files.
 *
 * Now supports additional operations to seek and query the next contig from the file, for efficient
 * implementation of jumping forward in the file.
 *
 * @author Tim Fennell
 * @author Extended in parts by Mark DePristo
 */
public class FastaSequenceFile2 implements ReferenceSequenceFile {
    private final File file;
    private BufferedInputStream in;
    private SAMSequenceDictionary sequenceDictionary = null;
    private String currentContigName = null;

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(FastaSequenceFile2.class);
    
    /**
     * Set to true to see lots of debugging output during operation
     */
    private final boolean DEBUG = false;

    /**
     * The name, if known, of the next contig in the file.  Can be null, indicating either that there is
     * no known next contig name, or that the file has been completely read
     */
    private String nextContigName = null;

    /** Constructs a FastaSequenceFile that reads from the specified file. */
    public FastaSequenceFile2(final File file) {
        this.file = file;
        initializeInputStream();

        // Try and locate the dictionary
        String dictionaryName = file.getAbsolutePath();
        dictionaryName = dictionaryName.substring(0, dictionaryName.lastIndexOf(".fasta"));
        dictionaryName += ".dict";
        final File dictionary = new File(dictionaryName);
        if (dictionary.exists()) {
            IoUtil.assertFileIsReadable(dictionary);

            try {
                final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                final SAMFileHeader header = codec.decode(new AsciiLineReader(new FileInputStream(dictionary)), dictionary);
                if (header.getSequenceDictionary() != null && header.getSequenceDictionary().size() > 0) {
                    this.sequenceDictionary = header.getSequenceDictionary();
                }
            }
            catch (Exception e) {
                throw new PicardException("Could not open sequence dictionary file: " + dictionaryName, e);
            }
        }
    }

    /**
     * Opens the file for input.  Closes the previous input stream if there is one, and sets up all of the key
     * variables such that this object's state is clean and ready for processing sequences.
     */
    private void initializeInputStream() {
        // we're at the start and haven't processed anything yet
        nextContigName = currentContigName = null;

        if ( this.in != null )  { // this isn't our first time here
            try {
                this.in.close();
            } catch ( IOException e ) {
                throw new RuntimeIOException("initializing InputStream failure", e);
            }
        }

        // Now reopen the input stream
        this.in = new BufferedInputStream(IoUtil.openFileForReading(file));
    }

    /**
     * Returns the list of sequence records associated with the reference sequence if found
     * otherwise null.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return this.sequenceDictionary;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Support functions for seeking around in the file
    //
    // --------------------------------------------------------------------------------------------------------------
    
    /**
     * Returns the distance, in bp, between contig1 and contig2 according to this fasta file's dictionary.  That is,
     * the number of bp we'd have to traverse to move from the start of contig1 to reach the start of contig2.
     *
     * If contig1 occurs before contig2, a negative number is returned.  0 indicates that the contigs are the same.
     *
     * Returns Integer.MAX_VALUE if the sequence dictionary cannot be found
     *
     * @param contig1
     * @param contig2
     * @return distance in bp from the start of contig1 to the start of contig2
     */
    public long getDistanceBetweenContigs(final String contig1, final String contig2) {
        assert contig1 != null;
        assert contig2 != null;

        final SAMSequenceDictionary seqDict = getSequenceDictionary();

        if ( seqDict == null ) // we couldn't load the reference dictionary
            return Integer.MAX_VALUE;

        SAMSequenceRecord contig1Rec = seqDict.getSequence(contig1);
        SAMSequenceRecord contig2Rec = seqDict.getSequence(contig2);

        assert contig1Rec != null : "Contig1 record is null: " + contig1;
        assert contig1Rec != null : "Contig2 record is null: " + contig2;

        if ( DEBUG )
            logger.debug(String.format("Contig1=(%s, %d), contig2=(%s, %d)%n",
                    contig1, contig1Rec.getSequenceIndex(),
                    contig2, contig2Rec.getSequenceIndex()));

        int startIndex = Math.min(contig1Rec.getSequenceIndex(), contig2Rec.getSequenceIndex());
        int lastIndex = Math.max(contig1Rec.getSequenceIndex(), contig2Rec.getSequenceIndex());

        long bytesToTraverse = 0;
        for ( int i = startIndex; i < lastIndex; i++ ) {
            SAMSequenceRecord rec = seqDict.getSequence(i);
            bytesToTraverse += rec.getSequenceLength();
            if ( DEBUG )
                logger.debug(String.format("  -> Traversing from %15s to %15s requires reading at least %10d bytes to pass contig %15s, total bytes %10d%n",
                    contig1, contig2, rec.getSequenceLength(), rec.getSequenceName(), bytesToTraverse));
        }

        if ( contig1Rec.getSequenceIndex() > contig2Rec.getSequenceIndex() )
            bytesToTraverse *= -1;  // we are going backward!

        if ( DEBUG ) logger.debug(String.format("  -> total distance is %d%n", bytesToTraverse));

        return bytesToTraverse;
    }

    /**
     * Seeks to seekContig in the fasta file, such that nextSequence() will read the seekContig from the fasta
     * file.  Only allows forward seeks.  Throws a RuntimeIOException if the seekContig is before the current
     * contig.
     *
     * @param seekContig the contig I want to seek to
     * @return true on success
     *
     * @see #seekToContig(String)
     */
    public boolean seekToContig( final String seekContig ) {
        return seekToContig(seekContig, false);
    }

    /**
     * Seeks to seekContig in the fasta file, such that nextSequence() will read the seekContig from the fasta
     * file.  If enableBacktracking is false, only allows forward seeks, and throws a RuntimeIOException if the
     * seekContig is before the current contig.  If enableBacktracking is true, then if seekContig is before
     * the current contig, resets the input stream and seeks to the contig.
     *
     * Requires that the fasta file have a SequenceDictionary associated with it.  Otherwises throws an error
     *
     * @param seekContig The contig I want to seek to
     * @param enableBacktracking Should we allow seeks to contigs earlier in the file?
     * @return true on success
     */
    public boolean seekToContig(final String seekContig, boolean enableBacktracking ) {
        if ( DEBUG ) logger.debug(String.format("seekToContig( %s, %b )%n", seekContig, enableBacktracking));

        String curContig = getContigName();
        String nextContig = null;
        
        if ( curContig == null ) {
            logger.info(String.format("CurrentContig is null"));
            if ( this.sequenceDictionary == null )
                throw new PicardException( String.format("Seeking within contigs requires FASTA dictionary, but none was available for %s", this.file ));

            // We only reach this point when we're seeking before we've read in any of the fasta file,
            // so assume we are at the start of the file
            nextContig = this.sequenceDictionary.getSequence(0).getSequenceName();
        }
        else
            nextContig = getNextContigName();

        if ( nextContig == null )   // we're are at the end of the stream
            return false;

        // we have already read in the current contig, we are jumping from the next contig onwards
        long dist = getDistanceBetweenContigs(nextContig, seekContig);

        if ( dist == Integer.MAX_VALUE )
            return false;       // we don't know where to go
        else if ( dist == 0 )
            return true;        // We already here!
        else if ( dist < 0 ) {
            if ( enableBacktracking ) {
                // System.out.printf("*** Backtracking to %s%n", seekContig);
                // restart from the beginning, and try again
                initializeInputStream();
                return seekToContig(seekContig, enableBacktracking);
            } else
                return false;       // we're not going backwards just yet
        }
        else {
            if ( DEBUG ) logger.debug(String.format("Going to seek to contig %s with skip %d%n", seekContig, dist));
            // we're actually going to jump somewhere, so prepare the state
            this.nextContigName = null;         // reset the contig info

            // TODO: this is a dangerous method -- can we get access to the underlying file object seek?
            long bytesToSkip = dist;
            while ( bytesToSkip > 0 ) {
                try {
                    final long skipped = this.in.skip(bytesToSkip);
                    bytesToSkip -= skipped;
                    // System.out.printf("  -> skipping %d, %d remaining%n", skipped, bytesToSkip);
                } catch (IOException ioe) {
                    throw new PicardException("Error reading from file: " + this.file.getAbsolutePath(), ioe);
                }
            }

            if ( bytesToSkip != 0 ) { // skip dist bytes
                throw new PicardException(String.format("Failed to skip all of the %d bases requested, only got %d", dist, dist - bytesToSkip * 2));
            }

            // at this point we're ready to start looking for the next header, so call seekNextContigName()
            final String next = seekForNextContig(seekContig);

            if ( ! next.equals(seekContig) ) // OMG, what the hell happened, throw a runtime exception
                throw new PicardException(String.format("Failed to seek from %s to %s, ended up at %s",
                            curContig, seekContig, next));
            else {
                this.currentContigName = next;
                return true;
            }
        }
    }

    /**
     * Reads the next contig from the fasta file, and returns it as a ReferenceSequence.
     *  
     * @return null if there are no more sequences in the fasta stream
     */
    public ReferenceSequence nextSequence() {
        if ( DEBUG ) logger.debug(String.format("Calling nextSequence()%n"));
        
        // Read the header line
        currentContigName = getNextContigName();
        if ( currentContigName == null ) return null; // no more sequences!

        int index = -1;
        if ( this.sequenceDictionary != null )
           index = this.sequenceDictionary.getSequenceIndex(currentContigName);

        // Read the sequence
        byte[] tmp = new byte[4096];
        int basesRead;
        int totalBasesRead = 0;
        final int knownLength = (index == -1) ? -1 : this.sequenceDictionary.getSequence(index).getSequenceLength();
        final int lengthByteArray = (knownLength != -1) ? knownLength : 250000000;
        byte[] bases = new byte[lengthByteArray];

        while ((basesRead = readNextLine(bases, totalBasesRead)) != 0) {
            totalBasesRead += basesRead;

            // Make sure we'll have space for the next iteration if we need it
            if (totalBasesRead == knownLength) {
                //System.out.printf("Read bases: %s%n", StringUtil.bytesToString(bases, totalBasesRead - basesRead, basesRead).trim());

                assert peekOneByte() == -1 || peekOneByte() == '>' : "We somehow managed to read in enough bytes for the contig, but didn't pass through the entire contig";
                break;
            } else {
                final byte b = peekOneByte();
                if (b == -1 || b == '>') {
                    break;
                }
                else if (totalBasesRead == bases.length) {
                    tmp = new byte[bases.length * 2];
                    System.arraycopy(bases, 0, tmp, 0, totalBasesRead);
                    bases = tmp;
                    tmp = null;
                }
            }
        }

        // And lastly resize the array down to the right size
        if (totalBasesRead != bases.length) {
            tmp = new byte[totalBasesRead];
            System.arraycopy(bases, 0, tmp, 0, totalBasesRead);
            bases = tmp;
            tmp = null;
        }

        assert knownLength == -1 || knownLength == bases.length;

        this.nextContigName = null; // we no longer know what the next contig name is

        if ( DEBUG ) logger.debug(String.format(" => nextSequence() is returning %s, known length = %d%n", this.currentContigName, knownLength));
        if ( DEBUG ) logger.debug(String.format(" => nextSequence() next is %s%n", this.getNextContigName()));

        return new ReferenceSequence(currentContigName, index, bases);
    }

    /**
     * Returns the next of the next contig, or null if there are no more contigs. Stateful function, it
     * remembers the name of the next contig.  Use readNextContigName for stateless operation.  This function
     * is primarily useful if you need to know what the next contig is in the stream for algorithmic purposes.
     *
     * Note that this function assumes the stream is sitting right at the end of the previous contig, or at the
     * beginning of the file.
     * 
     * @return the name of the next contig, or null if there is no next contig
     */
    public String getNextContigName() {
        if ( DEBUG ) logger.debug(String.format("getNextContigName() => %s%n", this.nextContigName));

        if ( this.nextContigName == null ) {
            // If it's not null, we've already looked up the next contig name, just return it and happily continue
            // Otherwise we need to actually read in the name
            this.nextContigName = readNextContigName();
        }

        if ( DEBUG ) logger.debug(String.format("nextContigName is now %s%n", nextContigName));
        return this.nextContigName;
    }

    /**
     * Simply reads the next contig name from the fasta input stream.  It assumes that the stream is positioned
     * immediately before the start of the next contig.  Calling it without this condition met results in the
     * RuntimeIOException being thrown.  This method advances the fasta stream itself -- subsequent calls to the
     * method will lead to errors.  getNextContigName is the stateful version.
     *
     * @See getNextContigName()
     *
     * @return the string name of the next contig
     */
    private String readNextContigName() {
        // Otherwise we need to actually read in the name
        byte[] tmp = new byte[4096];
        final int nameLength = readNextLine(tmp, 0);
        if (nameLength != 0) {
            // 0 means no more sequences!
            if ( tmp[0] != '>' )
                throw new RuntimeIOException("The next line is supposed to be a fasta contig start but found " + StringUtil.bytesToString(tmp, 0, nameLength).trim());
            
            return StringUtil.bytesToString(tmp, 1, nameLength).trim();
        }
        
        return null;
    }

    /**
     * @return The name of the contig we returned in the last call to nextSequence()
     */
    public String getContigName() {
        return this.currentContigName;
    }

    /**
     * Moves the IO stream to right before the next contig marker in the fasta file, for anywhere inside the
     * previous contig.  It is primarily useful as a supplementary routine for jumping forward in the file by
     * N bases, since we don't know exactly how far way N bases in bytes will be in the file.  So a guess jump
     * will put us somewhere before the target contig, and we use this routine to seek forward to the actual
     * contig we want.
     * 
     * @return the name of the next contig, as a string
     */
    public String seekForNextContig(final String targetContig ) {
        //System.out.printf("seekForNextContig()%n");

        int basesRead;
        int totalBasesRead = 0;
        byte[] bases = new byte[4096];
        int i = 0;
        while ((basesRead = readNextLine(bases, 0)) != 0) {
            totalBasesRead += basesRead;

            // Keep looking for the > marking the start of the line, and stop
            final byte b = peekOneByte();
            if (b == -1 || b == '>') {
                final String foundContig = readNextContigName();
                // System.out.printf("Found a contig name line %s%n", foundContig);
                final int foundIndex = this.sequenceDictionary.getSequenceIndex(foundContig);
                final int ourIndex = this.sequenceDictionary.getSequenceIndex(targetContig);

                if ( foundIndex == ourIndex ) {
                    // we found our target!
                    this.nextContigName = foundContig;   // store the right answer
                    if ( DEBUG ) logger.debug(String.format("seekForNextContig found %s%n", foundContig));
                    return foundContig;
                }
                else if ( foundIndex <= ourIndex )
                    // we are still looking for our contig, the seek estimate was inaccurate relative to the size of contings in this area
                    continue;
                else {
                    // This is really bad -- we are past our target
                    throw new RuntimeIOException(String.format("Seek pushes us past our target contig of %s, instead we found %s, which is after the target in the sequence dictions", targetContig, foundContig));
                }
            }
        }


        return null;
    }

    /** Peeks one non line-terminating byte from the file and returns it. */
    private byte peekOneByte() {
        try {
            this.in.mark(16);
            byte b = '\n';
            while (b == '\n' || b == '\r') {
                b = (byte) this.in.read();
            }

            this.in.reset();
            return b;
        }
        catch (IOException ioe) {
            throw new PicardException("Error reading from file: " + this.file.getAbsolutePath(), ioe);
        }
    }

    /**
     * Reads the next line from the file and puts the bytes into the buffer starting at
     * the provided index. Stops when the buffer is full, a line terminator is hit or
     * the end of the file is hit.
     */
    private int readNextLine(final byte[] buffer, final int index) {
        try {
            int next;
            int read = 0;
            while ((next =  this.in.read()) != -1 && index + read < buffer.length) {
                final byte b = (byte) next;
                if (b == '\r' || b == '\n') {
                    if (read != 0) return read;
                }
                else {
                    buffer[index + read++] = b;
                }
            }
            return read;
        }
        catch (IOException ioe) {
            throw new PicardException("Error reading line from file: " + this.file.getAbsolutePath(), ioe);
        }
    }

    /** Returns the full path to the reference file. */
    public String toString() {
        return this.file.getAbsolutePath();
    }
}
