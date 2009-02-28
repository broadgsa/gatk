package edu.mit.broad.picard.reference;

import java.io.File;

/**
 * Factory class for creating ReferenceSequenceFile instances for reading reference
 * sequences store in various formats.
 *
 * @author Tim Fennell
 */
public class ReferenceSequenceFileFactory {

    /**
     * Attempts to determine the type of the reference file and return an instance
     * of ReferenceSequenceFile that is appropriate to read it.
     *
     * @param file the reference sequence file on disk
     */
    public static ReferenceSequenceFile getReferenceSequenceFile(File file) {
        String name = file.getName();
        if (name.endsWith(".fasta") || name.endsWith("fasta.gz") || name.endsWith(".txt") || name.endsWith(".txt.gz")) {
            return new FastaSequenceFile(file);
        }
        else {
            throw new IllegalArgumentException("File is not a supported reference file type: " + file.getAbsolutePath());
        }
    }
}
