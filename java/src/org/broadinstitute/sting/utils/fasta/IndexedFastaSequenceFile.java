package org.broadinstitute.sting.utils.fasta;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMTextHeaderCodec;
import net.sf.samtools.util.AsciiLineReader;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Iterator;

/**
 * A fasta file driven by an index for fast, concurrent lookups.  Supports two interfaces:
 * the ReferenceSequenceFile for old-style, stateful lookups and a direct getter.
 */
public class IndexedFastaSequenceFile implements ReferenceSequenceFile {
    /**
     * Stores the main fasta file.
     */
    private final File file;

    /**
     * The interface facilitating direct access to the fasta.
     */
    private FileChannel channel;

    /**
     * A representation of the sequence dictionary, stored alongside the fasta in a .dict file.
     */
    private SAMSequenceDictionary sequenceDictionary = null;

    /**
     * A representation of the sequence index, stored alongside the fasta in a .fasta.fai file.
     */
    private FastaSequenceIndex index;

    /**
     * An iterator into the fasta index, for traversing iteratively across the fasta.
     */
    private Iterator<FastaSequenceIndexEntry> indexIterator;

    /**
     * Open the given indexed fasta sequence file.  Throw an exception if the file cannot be opened.
     * @param file The file to open.
     * @throws FileNotFoundException If the fasta or any of its supporting files cannot be found.
     */
    public IndexedFastaSequenceFile(File file) throws FileNotFoundException {
        this.file = file;
        FileInputStream in = new FileInputStream(file);
        channel = in.getChannel();

        loadDictionary(file);
        loadIndex(file);
        sanityCheckDictionaryAgainstIndex();
    }

    /**
     * Loads a dictionary, if available.
     * @param fastaFile File to check for a match.
     */
    private void loadDictionary( File fastaFile ) {
        // Try and locate the dictionary
        String dictionaryName = fastaFile.getAbsolutePath();
        dictionaryName = dictionaryName.substring(0, getFastaFileExtensionStart(dictionaryName));
        dictionaryName += ".dict";
        final File dictionary = new File(dictionaryName);
        if (!dictionary.exists())
            throw new PicardException("Unable to load .dict file.  Dictionary is required for the indexed fasta reader.");    

        IoUtil.assertFileIsReadable(dictionary);

        try {
            final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            final SAMFileHeader header = codec.decode(new AsciiLineReader(new FileInputStream(dictionary)),
                                                      dictionary.toString());
            if (header.getSequenceDictionary() != null && header.getSequenceDictionary().size() > 0) {
                this.sequenceDictionary = header.getSequenceDictionary();
            }
        }
        catch (Exception e) {
            throw new PicardException("Could not open sequence dictionary file: " + dictionaryName, e);
        }

    }

    /**
     * Gets the index of the first character in the fasta file's extension.
     * @param filename The filename of the fasta.  Must not be null, and must end with either '.fasta' or '.fa'.
     * @return The index of the start of the extension within the filename.  If neither '.fasta' nor '.fa' are
     *         present in the filename, a StingException will be thrown.
     */
    private int getFastaFileExtensionStart( String filename ) {
        if( filename.endsWith(".fasta") )
            return filename.lastIndexOf(".fasta");
        else if( filename.endsWith(".fa") )
            return filename.lastIndexOf(".fa");
        else
            throw new StingException("Invalid fasta filename; fasta filename must end with '.fasta' or '.fa'.");
    }

    /**
     * Loads the index for the fasta, if present.  Throws an exception if now present.
     * @param fastaFile FASTA file to load.
     * @throws FileNotFoundException if FASTA file cannot be found.
     */
    private void loadIndex( File fastaFile ) throws FileNotFoundException {
        File indexFile = new File(fastaFile.getAbsolutePath() + ".fai");
        if (!indexFile.exists())
            throw new PicardException(String.format("Unable to load fasta index file %s.  "+
                                                    "Please create it using 'samtools faidx'.",indexFile.getAbsolutePath()));
        index = new FastaSequenceIndex(indexFile);
        reset();
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

        FastaSequenceIndexEntry indexEntry = index.getIndexEntry(contig);

        if(stop > indexEntry.getSize())
            throw new PicardException("Query asks for data past end of contig");

        int length = (int)(stop - start + 1);

        byte[] target = new byte[length];
        ByteBuffer targetBuffer = ByteBuffer.wrap(target);

        final int basesPerLine = indexEntry.getBasesPerLine();
        final int bytesPerLine = indexEntry.getBytesPerLine();

        final long startOffset = ((start-1)/basesPerLine)*bytesPerLine + (start-1)%basesPerLine;
        final long stopOffset = ((stop-1)/basesPerLine)*bytesPerLine + (stop-1)%basesPerLine;
        final int size = (int)(stopOffset-startOffset)+1;

        ByteBuffer channelBuffer = ByteBuffer.allocate(size);
        try {
            channel.read(channelBuffer,indexEntry.getLocation()+startOffset);
        }
        catch(IOException ex) {
            throw new PicardException("Unable to map FASTA file into memory.");
        }

        channelBuffer.position(0);
        channelBuffer.limit(Math.min(basesPerLine-(int)startOffset%bytesPerLine,size));

        while( channelBuffer.hasRemaining() ) {
            targetBuffer.put(channelBuffer);

            channelBuffer.limit(Math.min(channelBuffer.limit()+bytesPerLine,size));
            channelBuffer.position(Math.min(channelBuffer.position()+bytesPerLine-basesPerLine,size));
        }

        return new ReferenceSequence( contig, sequenceDictionary.getSequenceIndex(contig), target );
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

    /**
     * Reset the iterator over the index.
     */
    @Override
    public void reset() {
        indexIterator = index.iterator();
    }

    /**
     * A simple toString implementation for debugging.
     * @return String representation of the file.
     */
    public String toString() {
        return this.file.getAbsolutePath();
    }
}
