package org.broadinstitute.sting.utils.fasta;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;

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
import java.util.Iterator;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMTextHeaderCodec;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.AsciiLineReader;

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

    private SAMSequenceDictionary sequenceDictionary = null;    

    private FastaSequenceIndex index;
    private Iterator<FastaSequenceIndexEntry> indexIterator;

    public IndexedFastaSequenceFile(File file) throws FileNotFoundException {
        this.file = file;
        // TODO: Add support for gzipped files
        in = new FileInputStream(file);
        channel = in.getChannel();

        loadDictionary(file);
        loadIndex(file);
        sanityCheckDictionaryAgainstIndex();
    }

    /**
     * Loads a dictionary, if available.
     * @param fastaFile File to check for a match.
     * TODO: This code is copied directly from FastaSequenceFile / FastaSequenceFile2.  Bring it into a shared utility.
     */
    private void loadDictionary( File fastaFile ) {
        // Try and locate the dictionary
        String dictionaryName = fastaFile.getAbsolutePath();
        dictionaryName = dictionaryName.substring(0, dictionaryName.lastIndexOf(".fasta"));
        dictionaryName += ".dict";
        final File dictionary = new File(dictionaryName);
        if (!dictionary.exists())
            throw new PicardException("Unable to load .dict file.  Dictionary is required for the indexed fasta reader.");    

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

    /**
     * Loads the index for the fasta, if present.  Throws an exception if now present.
     */
    private void loadIndex( File fastaFile ) throws FileNotFoundException {
        File indexFile = new File(fastaFile.getAbsolutePath() + ".fai");
        if (!indexFile.exists())
            throw new PicardException(String.format("Unable to load fasta index file %s.  "+
                                                    "Please create it using 'samtools faidx'.",indexFile.getAbsolutePath()));
        index = new FastaSequenceIndex(indexFile);
        indexIterator = index.iterator();
    }

    /**
     * Do some basic checking to make sure the dictionary and the index match.
     */
    private void sanityCheckDictionaryAgainstIndex() {
        // Make sure dictionary and index are the same size.
        if( sequenceDictionary.getSequences().size() != index.size() )
            throw new PicardException("Sequence dictionary and index contain different numbers of contigs");

        for( SAMSequenceRecord sequenceRecord: sequenceDictionary.getSequences() ) {
            // Make sure sequence name is present in the index.
            String sequenceName = sequenceRecord.getSequenceName();
            if( !index.hasIndexEntry(sequenceName) )
                throw new PicardException("Index does not contain dictionary entry: " + sequenceName );

            // Make sure sequence length matches index length.
            if( sequenceRecord.getSequenceLength() != index.getIndexEntry(sequenceName).getSize())
                throw new PicardException("Index length does not match dictionary length for contig: " + sequenceName );
        }
    }

    /**
     * Retrieves the sequence dictionary for the fasta file.
     * @return sequence dictionary of the fasta.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }

    /**
     * Retrieves the complete sequence described by this contig.
     * @param contig contig whose data should be returned.
     * @return The full sequence associated with this contig.
     */
    public ReferenceSequence getSequence( String contig ) {
        return getSubsequenceAt( contig, 1, (int)index.getIndexEntry(contig).getSize() );
    }

    /**
     * Gets the subsequence of the contig in the range [start,stop]
     * @param contig Contig whose subsequence to retrieve.
     * @param start inclusive, 1-based start of region.
     * @param stop inclusive, 1-based stop of region.
     * @return The partial reference sequence associated with this range.
     */
    public ReferenceSequence getSubsequenceAt( String contig, long start, long stop ) {
        if(start > stop)
            throw new PicardException(String.format("Malformed query; start point %d lies after end point %d",start,stop));
        if(start > Integer.MAX_VALUE)
            throw new PicardException("Due to current ReferenceSequence limitations, a start point larger than Integer.MAX_VALUE cannot be loaded.");
        if(stop - start + 1 > Integer.MAX_VALUE)
            throw new PicardException("Due to current ReferenceSequence limitations, a region larger than Integer.MAX_VALUE cannot be loaded.");

        FastaSequenceIndexEntry indexEntry = index.getIndexEntry(contig);

        if(stop > indexEntry.getSize())
            throw new PicardException("Query asks for data past end of contig");

        int length = (int)(stop - start + 1);

        final int basesPerLine = indexEntry.getBasesPerLine();
        final int bytesPerLine = indexEntry.getBytesPerLine();

        // Start reading at the closest start-of-line to our data.
        long readStart = indexEntry.getLocation() + ((start-1) / basesPerLine) * bytesPerLine;
        int dataOfInterestStart = (int)((start-1) % basesPerLine);

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

        return new ReferenceSequence( contig, sequenceDictionary.getSequenceIndex(contig), accumulator );
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

    /**
     * Gets the next sequence if available, or null if not present.
     * @return next sequence if available, or null if not present.
     */
    public ReferenceSequence nextSequence() {
        if( !indexIterator.hasNext() )
            return null;
        return getSequence( indexIterator.next().getContig() );
    }

    public String toString() {
        return this.file.getAbsolutePath();
    }
}
