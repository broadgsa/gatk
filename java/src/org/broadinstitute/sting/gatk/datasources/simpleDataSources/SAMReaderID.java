package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import java.io.File;
import java.util.List;
import java.util.Collections;

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
     * A list of tags associated with this BAM file.
     */
    protected final List<String> tags;

    /**
     * Creates an identifier for a SAM file based on read.
     * @param samFile The source file for SAM data.
     * @param tags tags to use when creating a reader ID.
     */
    public SAMReaderID(File samFile, List<String> tags) {
        this.samFile = samFile;
        this.tags = tags;
    }

    /**
     * Creates an identifier for a SAM file based on read.
     * @param samFileName The source filename for SAM data.
     * @param tags tags to use when creating a reader ID.
     */
    public SAMReaderID(String samFileName, List<String> tags) {
        this(new File(samFileName),tags);        
    }

    /**
     * Gets the tags associated with the given BAM file.
     * @return A collection of the tags associated with this file.
     */
    public List<String> getTags() {
        return Collections.unmodifiableList(tags);
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
