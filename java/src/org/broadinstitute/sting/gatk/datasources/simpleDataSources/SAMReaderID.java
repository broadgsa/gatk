package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import java.io.File;

/**
 * Uniquely identifies a SAM file reader.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMReaderID {
    /**
     * The SAM file at the heart of this reader.  SAMReaderID
     * currently supports only file-based readers.
     */
    protected final File samFile;

    /**
     * Creates an identifier for a SAM file based on read.
     * @param samFile The source file for SAM data.
     */
    protected SAMReaderID(File samFile) {
        this.samFile = samFile;
    }

    /**
     * Compare two IDs to see whether they're equal.
     * @param other The other identifier.
     * @return True iff the two readers point to the same file.
     */
    public boolean equals(Object other) {
        if(other == null) return false;
        if(!(other instanceof SAMReaderID)) return false;

        SAMReaderID otherID = (SAMReaderID)other;
        return this.samFile.equals(otherID.samFile);
    }

    /**
     * Generate a hash code for this object.
     * @return A hash code, based solely on the file name at this point.
     */
    public int hashCode() {
        return samFile.hashCode();
    }
}
