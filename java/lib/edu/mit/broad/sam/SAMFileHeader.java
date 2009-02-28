/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.sam;


import java.util.*;

/**
 * Header information from a SAM file.
 */
public class SAMFileHeader
{
    public static final String VERSION_TAG = "VN";
    public static final String CURRENT_VERSION = "1.0";

    public enum SortOrder {

        unsorted(null),
        queryname(SAMRecordQueryNameComparator.class),
        coordinate(SAMRecordCoordinateComparator.class);

        private Class<? extends SAMRecordComparator> comparator;

        SortOrder(final Class<? extends SAMRecordComparator> comparatorClass) {
            this.comparator = comparatorClass;
        }

        public Class<? extends SAMRecordComparator> getComparator() {
            return comparator;
        }
    }

    public enum GroupOrder {
        none, query, reference
    }

    private final Map<String, Object> mAttributes =
        new HashMap<String, Object>();
    private List<SAMSequenceRecord> mSequences =
        new ArrayList<SAMSequenceRecord>();
    private List<SAMReadGroupRecord> mReadGroups =
        new ArrayList<SAMReadGroupRecord>();
    private final List<SAMProgramRecord> mProgramRecords = new ArrayList<SAMProgramRecord>();
    private final Map<String, SAMSequenceRecord> mSequenceMap =
        new HashMap<String, SAMSequenceRecord>();
    private final Map<String, SAMReadGroupRecord> mReadGroupMap =
        new HashMap<String, SAMReadGroupRecord>();
    private Map<String, SAMProgramRecord> mProgramRecordMap = new HashMap<String, SAMProgramRecord>();

    public SAMFileHeader() {
        setAttribute(VERSION_TAG, CURRENT_VERSION);
    }

    public String getVersion() {
        return (String) getAttribute("VN");
    }

    public String getCreator() {
        return (String) getAttribute("CR");
    }

    public Object getAttribute(final String key) {
        return mAttributes.get(key);
    }

    public Set<Map.Entry<String, Object>> getAttributes() {
        return mAttributes.entrySet();
    }

    public List<SAMSequenceRecord> getSequences() {
        return mSequences;
    }

    public List<SAMReadGroupRecord> getReadGroups() {
        return mReadGroups;
    }

    public SAMSequenceRecord getSequence(final String name) {
        return mSequenceMap.get(name);
    }

    public SAMReadGroupRecord getReadGroup(final String name) {
        return mReadGroupMap.get(name);
    }

    public void setSequences(final List<SAMSequenceRecord> list) {
        mSequences = list;
        mSequenceMap.clear();
        int index = 0;
        for (final SAMSequenceRecord record : list) {
            record.setSequenceIndex(index++);
            mSequenceMap.put(record.getSequenceName(), record);
        }
    }

    public SAMSequenceRecord getSequence(final int sequenceIndex) {
        if (sequenceIndex < 0 || sequenceIndex >= mSequences.size()) {
            return null;
        }
        return mSequences.get(sequenceIndex);
    }

    public int getSequenceIndex(final String sequenceName) {
        final SAMSequenceRecord record = mSequenceMap.get(sequenceName);
        if (record == null) {
            return -1;
        }
        return record.getSequenceIndex();
    }

    public void setAttribute(final String key, final String value) {
        mAttributes.put(key, value);
    }

    public void setReadGroups(final List<SAMReadGroupRecord> readGroups) {
        mReadGroups = readGroups;
        mReadGroupMap.clear();
        for (final SAMReadGroupRecord readGroupRecord : readGroups) {
            mReadGroupMap.put(readGroupRecord.getReadGroupId(), readGroupRecord);
        }
    }

    public List<SAMProgramRecord> getProgramRecords() {
        return Collections.unmodifiableList(mProgramRecords);
    }

    public void addProgramRecord(final SAMProgramRecord programRecord) {
        this.mProgramRecords.add(programRecord);
        this.mProgramRecordMap.put(programRecord.getProgramGroupId(), programRecord);
    }

    public SAMProgramRecord getProgramRecord(final String name) {
        return this.mProgramRecordMap.get(name);
    }

    public SortOrder getSortOrder() {
        if (getAttribute("SO") == null) {
            return SortOrder.unsorted;
        }
        return SortOrder.valueOf((String)getAttribute("SO"));
    }

    public void setSortOrder(final SortOrder so) {
        setAttribute("SO", so.name());
    }

    public GroupOrder getGroupOrder() {
        if (getAttribute("GO") == null) {
            return GroupOrder.none;
        }
        return GroupOrder.valueOf((String)getAttribute("GO"));
    }

    public void setGroupOrder(final GroupOrder go) {
        setAttribute("GO", go.name());
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SAMFileHeader that = (SAMFileHeader) o;

        if (mAttributes != null ? !mAttributes.equals(that.mAttributes) : that.mAttributes != null) return false;
        if (mProgramRecords != null ? !mProgramRecords.equals(that.mProgramRecords) : that.mProgramRecords != null)
            return false;
        if (mReadGroups != null ? !mReadGroups.equals(that.mReadGroups) : that.mReadGroups != null) return false;
        if (mSequences != null ? !mSequences.equals(that.mSequences) : that.mSequences != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = mAttributes != null ? mAttributes.hashCode() : 0;
        result = 31 * result + (mSequences != null ? mSequences.hashCode() : 0);
        result = 31 * result + (mReadGroups != null ? mReadGroups.hashCode() : 0);
        result = 31 * result + (mReadGroupMap != null ? mReadGroupMap.hashCode() : 0);
        result = 31 * result + (mProgramRecords != null ? mProgramRecords.hashCode() : 0);
        return result;
    }
}
