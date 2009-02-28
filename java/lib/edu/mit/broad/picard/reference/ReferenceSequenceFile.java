package edu.mit.broad.picard.reference;

import edu.mit.broad.sam.SAMSequenceRecord;

import java.util.List;

/**
 * An interface for working with files of reference sequences regardless of the file format
 * being used.
 *
 * @author Tim Fennell
 */
public interface ReferenceSequenceFile {

    /**
     * Must return a sequence dictionary with at least the following fields completed
     * for each sequence: name, length.
     *
     * @return a list of sequence records representing the sequences in this reference file
     */
    public List<SAMSequenceRecord> getSequenceDictionary();

    /**
     * Retrieves the next whole sequences from the file.
     * @return a ReferenceSequence or null if at the end of the file
     */
    public ReferenceSequence nextSequence();

}
