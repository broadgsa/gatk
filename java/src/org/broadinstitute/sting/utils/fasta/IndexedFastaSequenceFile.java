package org.broadinstitute.sting.utils.fasta;

import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.PicardException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.nio.ByteBuffer;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.Charset;
import java.nio.charset.CharacterCodingException;
import java.util.Scanner;

import net.sf.samtools.SAMSequenceDictionary;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 14, 2009
 * Time: 2:14:26 PM
 *
 * A fasta file driven by an index for fast, concurrent lookups.  Supports two interfaces:
 * the ReferenceSequenceFile for old-style, stateful lookups and a direct getter.
 */
public class IndexedFastaSequenceFile implements ReferenceSequenceFile {
    // Using buffer size of 4k because that's what Picard uses; no thought went into this.
    private static final int BUFFERSIZE = 4096;

    private final File file;
    private FileInputStream in;
    private FileChannel channel;

    private final FastaSequenceIndex index;

    private String currentContigName = null;

    public IndexedFastaSequenceFile(File file) throws FileNotFoundException {
        this.file = file;
        // TODO: Add support for gzipped files
        in = new FileInputStream(file);
        channel = in.getChannel();

        File indexFile = new File(file.getAbsolutePath() + ".fai");
        index = new FastaSequenceIndex(indexFile);
    }

    public SAMSequenceDictionary getSequenceDictionary() {
        throw new UnsupportedOperationException("Indexed fasta files do not require dictionaries");
    }

    public ReferenceSequence getSequence( String contig ) {
        return getSubsequenceAt( contig, 0, (int)index.getIndexEntry(contig).getSize() );
    }

    public ReferenceSequence getSubsequenceAt( String contig, int pos, int length ) {
        FastaSequenceIndexEntry indexEntry = index.getIndexEntry(contig);

        final int basesPerLine = indexEntry.getBasesPerLine();

        // Start reading at the closest start-of-line to our data.
        long readStart = indexEntry.getLocation() + (pos / basesPerLine);
        int dataOfInterestStart = pos % basesPerLine;

        byte[] accumulator = new byte[length];
        int nextAccumulatorSlot = 0;        

        while(length > 0) {
            ByteBuffer buffer = ByteBuffer.allocateDirect(BUFFERSIZE);
            try {
                channel.read(buffer, readStart);
                readStart += BUFFERSIZE;
            }
            catch( IOException ex ) {
                throw new PicardException("Unable to read directly from fasta", ex);
            }

            final int basesTransferred = transferToBuffer( buffer,
                                                           dataOfInterestStart,
                                                           accumulator,
                                                           nextAccumulatorSlot,
                                                           length );

            nextAccumulatorSlot += basesTransferred;
            length -= basesTransferred;
            dataOfInterestStart = 0;
        }

        return new ReferenceSequence( contig, pos, accumulator );
    }

    /**
     * Transfers the contents of the given ByteBuffer to the given byte array, discarding
     * line breaks at regular intervals.  Copies as many as length bases, depending on the
     * buffer size.  Returns the number of bytes actually copied.
     * @param source The source ByteBuffer.
     * @param sourceStart The starting position to copy within the byte buffer
     * @param target Destination for the data
     * @param targetStart Index into target buffer.
     * @param length How much data to move.
     * @return How many bytes were actually transferred.
     */
    private int transferToBuffer( ByteBuffer source,
                                  int sourceStart,
                                  byte[] target,
                                  int targetStart,
                                  int length ) {
        source.position(sourceStart);
        int basesRead = 0;
        CharsetDecoder decoder = Charset.forName("US-ASCII").newDecoder();

        Scanner scanner = null;
        try {
            scanner = new Scanner(decoder.decode(source).toString());
        }
        catch(CharacterCodingException ex) {
            throw new PicardException("Malformed subsequence",ex);
        }

        while( scanner.hasNext() && basesRead < length ) {
            String sourceLine = scanner.nextLine();
            byte[] sourceData = sourceLine.getBytes();
            int basesToTransfer = Math.min(sourceData.length,length - basesRead);
            System.arraycopy(sourceData,0,target,targetStart+basesRead,basesToTransfer);

            basesRead += basesToTransfer;
        }
        return basesRead;
    }


    public ReferenceSequence nextSequence() {
        return getSubsequenceAt("chrM", 0, 20);
    }

    public String toString() {
        return this.file.getAbsolutePath();
    }
}
