/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever.
* Neither the Broad Institute nor MIT can be responsible for its use, misuse, or
* functionality.
*/
package edu.mit.broad.picard.sam;

import edu.mit.broad.sam.*;
import edu.mit.broad.picard.PicardException;

import java.util.*;

/**
 * Merges SAMFileHeaders that have the same sequences into a single merged header
 * object while providing read group translation for cases where read groups
 * clash across input headers.
 *
 * @author Dave Tefft
 */
public class SamFileHeaderMerger {
    //Super Header to construct
    private final SAMFileHeader mergedHeader;
    private final Collection<SAMFileReader> readers;

    //Translation of old group ids to new group ids
    private final Map<SAMFileReader, Map<String, String>> samGroupIdTranslation =
            new HashMap<SAMFileReader, Map<String, String>>();

    //the groups from different files use the same group ids
    private boolean hasGroupIdDuplicates = false;

    //Translation of old program group ids to new program group ids
    private final Map<SAMFileReader, Map<String, String>> samProgramGroupIdTranslation =
            new HashMap<SAMFileReader, Map<String, String>>();

    //Letters to construct new ids from a counter
    private static final String ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";


    /**
     * Create SAMFileHeader with additional information
     *
     * @param readers same file readers to combine
     * @param sortOrder sort order new header should have
     */
    public SamFileHeaderMerger(final Collection<SAMFileReader> readers, final SAMFileHeader.SortOrder sortOrder) {
        this.readers = readers;
        this.mergedHeader = new SAMFileHeader();

        // Set sequences first because if it throws exception there is no need to continue
        final List<SAMSequenceRecord> sequences = getSAMSequences(readers);
        this.mergedHeader.setSequences(sequences);

        // Set program that creates input alignments
        for (final SAMProgramRecord program : mergeSAMProgramRecordLists(readers)) {
            this.mergedHeader.addProgramRecord(program);
        }

        // Set read groups for merged header
        final List<SAMReadGroupRecord> readGroups = getReadGroups(readers);
        this.mergedHeader.setReadGroups(readGroups);
        this.mergedHeader.setGroupOrder(SAMFileHeader.GroupOrder.none);

        this.mergedHeader.setSortOrder(sortOrder);
    }

    /**
     * Checks to see if there are clashes where different readers are using the same read
     * group IDs. If they are then a new set of unique read group IDs are generated (across all
     * read groups) otherwise the original read group headers are returned.
     *
     * @param readers readers to combine
     * @return new list of readgroups constructed from all the readers
     */
    private List<SAMReadGroupRecord> getReadGroups(final Collection<SAMFileReader> readers) {
        // Read groups as read from the readers
        final List<SAMReadGroupRecord> orginalReadGroups = new ArrayList<SAMReadGroupRecord>();

        // Read group with new ids that don't confict
        final List<SAMReadGroupRecord> modifiedReadGroups = new ArrayList<SAMReadGroupRecord>();

        //set to see if there are duplicate group ids and whether or not we need to modify them
        final Set<String> groupIdsSeenBefore = new HashSet<String>();

        int x = 0;
        this.hasGroupIdDuplicates = false;

        for (final SAMFileReader reader : readers) {
            final SAMFileHeader header = reader.getFileHeader();
            final Map<String, String> idTranslation = new HashMap<String, String>();

            // Iterate over read groups to find conflicting ids
            for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
                final String groupId = readGroup.getReadGroupId();
                final String newGroupId = createNewId(x++);

                // Check to see if same group id is used in two different readers
                if (groupIdsSeenBefore.contains(groupId)) {
                    hasGroupIdDuplicates = true;
                }
                groupIdsSeenBefore.add(groupId);

                // Creates a new read group with the new id and copies all it's attributes
                final SAMReadGroupRecord groupRecordWithNewId = copyReadGroup(readGroup, newGroupId);

                orginalReadGroups.add(readGroup);
                modifiedReadGroups.add(groupRecordWithNewId);

                idTranslation.put(groupId, newGroupId);
            }

            // Add id tranlation for updating SamRecords with new ids if neccessary
            this.samGroupIdTranslation.put(reader, idTranslation);
        }

        // return approriate readgroups whether or not the new ids have to be used
        if (this.hasGroupIdDuplicates) {
            return modifiedReadGroups;
        }
        else {
            return orginalReadGroups;
        }
    }

    /**
     * Get the sequences off the SAMFileReader header.  Throws runtime exception if the sequence
     * are different from one another
     *
     * @param readers readers to pull sequences from
     * @return sequences from files.  Each file should have the same sequence
     */
    private List<SAMSequenceRecord> getSAMSequences(final Collection<SAMFileReader> readers) {
        List<SAMSequenceRecord> sequences = null;
        for (final SAMFileReader reader : readers) {
            final SAMFileHeader header = reader.getFileHeader();

            if (sequences == null) {
                sequences = header.getSequences();
            }
            else {
                final List<SAMSequenceRecord> currentSequences = header.getSequences();
                if (!sequenceListsEqual(sequences, currentSequences)) {
                    throw new PicardException("Files are not compatible with each other.  They can not be combined");
                }
            }
        }
        return sequences;
    }

    /**
     * Checks the equality of two lists of sequence records using the isSameSequence
     * method instead of the equals method which is a more strict identity check.
     * @param s1 a list of sequence headers
     * @param s2 a second list of sequence headers
     * @return true if the two lists match otherwise false
     */
    private boolean sequenceListsEqual(final List<SAMSequenceRecord> s1, final List<SAMSequenceRecord> s2) {
        if (s1.size() != s2.size()) {
            return false;
        }
        for (int i = 0; i < s1.size(); ++i) {
            if (!s1.get(i).isSameSequence(s2.get(i))) {
                return false;
            }
        }
        return true;
    }

    /**
     * Find the alignment program that produced the readers.  If there are more than one
     * generate a new program represents that
     *
     * @param readers SAMFileReaders to pull program information from
     * @return SAMProgram record that represents all the readers
     */
    // TODO: this needs to be fixed up to support multiple program records (PIC-15)
    private List<SAMProgramRecord> mergeSAMProgramRecordLists(final Collection<SAMFileReader> readers) {
        final boolean programMixed = false;
        final List<SAMProgramRecord> ret = new ArrayList<SAMProgramRecord>();
        int nextProgramGroupId = 0;
        for (final SAMFileReader reader : readers) {
            final SAMFileHeader header = reader.getFileHeader();
            final Map<String, String> idTranslation = new HashMap<String, String>();
            for (final SAMProgramRecord oldProgramRecord : header.getProgramRecords()) {
                boolean foundMatch = false;
                for (final SAMProgramRecord newProgramRecord : ret) {
                    if (newProgramRecord.equivalent(oldProgramRecord)) {
                        idTranslation.put(oldProgramRecord.getProgramGroupId(), newProgramRecord.getProgramGroupId());
                        foundMatch = true;
                        break;
                    }
                }
                if (!foundMatch) {
                    final SAMProgramRecord newProgramRecord = new SAMProgramRecord(Integer.toString(nextProgramGroupId++));
                    copyProgramGroupAttributes(oldProgramRecord, newProgramRecord);
                    ret.add(newProgramRecord);
                    idTranslation.put(oldProgramRecord.getProgramGroupId(), newProgramRecord.getProgramGroupId());
                }
            }
            samProgramGroupIdTranslation.put(reader, idTranslation);
        }
        return ret;
    }

    private void copyProgramGroupAttributes(final SAMProgramRecord oldProgramRecord, final SAMProgramRecord newProgramRecord) {
        for (final Map.Entry<String, String> entry : oldProgramRecord.getAttributes()) {
            newProgramRecord.setAttribute(entry.getKey(), entry.getValue());
        }
    }


    /**
     * Copies all the attribute of a readgroup to a new readgroup with a new id
     *
     * @param readGroup  the group to be copied
     * @param modifiedId the id for the new readgroup
     * @return new read group
     */
    private SAMReadGroupRecord copyReadGroup(final SAMReadGroupRecord readGroup, final String modifiedId) {
        final SAMReadGroupRecord retval = new SAMReadGroupRecord(modifiedId);
        retval.setLibrary(readGroup.getLibrary());
        retval.setSample(readGroup.getSample());

        for (final Map.Entry<String, Object> attr : readGroup.getAttributes()) {
            retval.setAttribute(attr.getKey(), attr.getValue());
        }

        return retval;
    }


    /**
     * Creates a base 26 representation of an int
     *
     * @param n int to covert to letter representation
     * @return string rep for an int eg 0 = A  27 = AB
     */
    protected static String createNewId(int n) {
        final int base = ALPHABET.length();

        String s = "";
        while (true) {
            final int r = n % base;
            s = ALPHABET.charAt(r) + s;
            n = n / base;
            if (n == 0) {
                return s;
            }
            n -= 1;
        }
    }

    /** Returns the read group id that should be used for the input read and RG id. */
    public String getReadGroupId(final SAMFileReader reader, final String originalReadGroupId) {
        return this.samGroupIdTranslation.get(reader).get(originalReadGroupId);
    }

    /**
     * @param reader one of the input files
     * @param originalProgramGroupId a program group ID from the above input file
     * @return new ID from the merged list of program groups in the output file
     */
    public String getProgramGroupId(final SAMFileReader reader, final String originalProgramGroupId) {
        return this.samProgramGroupIdTranslation.get(reader).get(originalProgramGroupId);
    }

    /** Returns true if there are read group duplicates within the merged headers. */
    public boolean hasGroupIdDuplicates() {
        return this.hasGroupIdDuplicates;
    }

    /** Returns the merged header that should be written to any output merged file. */
    public SAMFileHeader getMergedHeader() {
        return this.mergedHeader;
    }

    /** Returns the collection of readers that this header merger is working with. */
    public Collection<SAMFileReader> getReaders() {
        return this.readers;
    }
}
