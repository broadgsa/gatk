package org.broadinstitute.sting.utils.fasta;

import edu.mit.broad.picard.PicardException;

import java.util.Scanner;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.MatchResult;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.File;
import java.io.FileNotFoundException;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 14, 2009
 * Time: 10:02:10 AM
 *
 * Reads a fasta index file (.fai).
 */
public class FastaSequenceIndex {
    private Map<String,FastaSequenceIndexEntry> sequenceEntries = 
            new HashMap<String,FastaSequenceIndexEntry>();

    /**
     * Build a sequence index from the specified file.
     * @param indexFile File to open.
     * @throws PicardException if file is of invalid format.
     */
    public FastaSequenceIndex( File indexFile ) throws FileNotFoundException {
        Scanner scanner = new Scanner(indexFile);

        while( scanner.hasNext() ) {
            // Tokenize and validate the index line.
            String result = scanner.findInLine("(\\w+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)");
            if( result == null )
                throw new PicardException("Found invalid line in index file:" + scanner.nextLine());
            MatchResult tokens = scanner.match();
            if( tokens.groupCount() != 5 )
                throw new PicardException("Found invalid line in index file:" + scanner.nextLine());

            // Skip past the line separator
            scanner.nextLine();

            // Parse the index line.
            String contig = tokens.group(1);
            long size = Long.valueOf(tokens.group(2));
            long location = Long.valueOf(tokens.group(3));
            int basesPerLine = Integer.valueOf(tokens.group(4));
            int bytesPerLine = Integer.valueOf(tokens.group(5));

            // Build sequence structure
            sequenceEntries.put( contig,new FastaSequenceIndexEntry(contig,location,size,basesPerLine,bytesPerLine) );
        }
    }

    /**
     * Does the given contig name have a corresponding entry?
     * @param contigName The contig name for which to search.
     * @return True if contig name is present; false otherwise.
     */
    public boolean hasIndexEntry( String contigName ) {
        return sequenceEntries.containsKey(contigName);
    }

    /**
     * Retrieve the index entry associated with the given contig.
     * @param contigName Name of the contig for which to search.
     * @return Index entry associated with the given contig.
     * @throws PicardException if the associated index entry can't be found.
     */
    public FastaSequenceIndexEntry getIndexEntry( String contigName ) {
        if( !hasIndexEntry(contigName) )
            throw new PicardException("Unable to find entry for contig: " + contigName);

        return sequenceEntries.get(contigName);
    }
}

class FastaSequenceIndexEntry {
    private String contig;
    private long location;
    private long size;
    private int basesPerLine;
    private int bytesPerLine;

    public FastaSequenceIndexEntry( String contig,
                                   long location,
                                   long size,
                                   int basesPerLine,
                                   int bytesPerLine ) {
        this.contig = contig;
        this.location = location;
        this.size = size;
        this.basesPerLine = basesPerLine;
        this.bytesPerLine = bytesPerLine;
    }

    /**
     * Gets the contig associated with this entry.
     * @return String representation of the contig.
     */
    public String getContig() {
        return contig;
    }

    /**
     * Gets the location of this contig within the fasta.
     * @return seek position within the fasta.
     */
    public long getLocation() {
        return location;
    }

    /**
     * Gets the size, in bytes, of the data in the contig.
     * @return size of the contig bases in bytes.
     */
    public long getSize() {
        return size;
    }

    /**
     * Gets the number of bases in a given line.
     * @return Number of bases in the fasta line.
     */
    public int getBasesPerLine() {
        return basesPerLine;
    }

    /**
     * How many bytes (bases + whitespace) are consumed by the
     * given line?
     * @return Number of bytes in a line.
     */
    public int getBytesPerLine() {
        return bytesPerLine;
    }

    public String toString() {
        return String.format("contig %s; location %d; size %d; basesPerLine %d; bytesPerLine %d", contig,
                                                                                                  location,
                                                                                                  size,
                                                                                                  basesPerLine,
                                                                                                  bytesPerLine );
    }
}
