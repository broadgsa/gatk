package org.broadinstitute.sting.utils.sam;

import java.util.List;

import net.sf.samtools.*;

/**
 * @author ebanks
 * GATKSAMRecord
 *
 * this class extends the samtools SAMRecord class and caches important
 * (and oft-accessed) data that's not already cached by the SAMRecord class
 *
 * IMPORTANT NOTE: Because ReadGroups are not set through the SAMRecord,
 *   if they are ever modified externally then one must also invoke the
 *   setReadGroup() method here to ensure that the cache is kept up-to-date.
 *
 */
public class GATKSAMRecord extends SAMRecord {

    // the underlying SAMRecord which we are wrapping
    private SAMRecord mRecord;

    // the SAMRecord data we're caching
    private String mReadString = null;
    private SAMReadGroupRecord mReadGroup = null;
    private Boolean mNegativeStrandFlag = null;
    private Boolean mUnmappedFlag = null;

    // because some values can be null, we don't want to duplicate effort
    private boolean retrievedReadGroup = false;

    public GATKSAMRecord(SAMRecord record) {
        super(null); // it doesn't matter - this isn't used
        if ( record == null )
            throw new IllegalArgumentException("The SAMRecord argument cannot be null");
        mRecord = record;

        // because attribute methods are declared to be final (and we can't overload them),
        // we need to actually set all of the attributes here

        // TODO -- remove this exception catch when Picard fixes its internal bug (PIC-291)

        try {
            List<SAMTagAndValue> attributes = record.getAttributes();
            for ( SAMTagAndValue attribute : attributes )
                setAttribute(attribute.tag, attribute.value);
        } catch (NullPointerException picardError) {
            ; // do nothing
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    // *** The following methods are overloaded to cache the appropriate data ***//
    ///////////////////////////////////////////////////////////////////////////////

    public String getReadString() {
        if ( mReadString == null )
            mReadString = mRecord.getReadString();
        return mReadString;
    }

    public void setReadString(String s) {
        mRecord.setReadString(s);
        mReadString = s;
    }

    public SAMReadGroupRecord getReadGroup() {
        if ( !retrievedReadGroup ) {
            SAMReadGroupRecord tempReadGroup = mRecord.getReadGroup();
            mReadGroup = (tempReadGroup == null ? tempReadGroup : new GATKSAMReadGroupRecord(tempReadGroup));
            retrievedReadGroup = true;
        }
        return mReadGroup;
    }

    public void setReadGroup(SAMReadGroupRecord record) {
        mReadGroup = record;
    }

    public boolean getReadUnmappedFlag() {
        if ( mUnmappedFlag == null )
            mUnmappedFlag = mRecord.getReadUnmappedFlag();
        return mUnmappedFlag;
    }

    public void setReadUmappedFlag(boolean b) {
        mRecord.setReadUmappedFlag(b);
        mUnmappedFlag = b;
    }

    public boolean getReadNegativeStrandFlag() {
        if ( mNegativeStrandFlag == null )
            mNegativeStrandFlag = mRecord.getReadNegativeStrandFlag();
        return mNegativeStrandFlag;
    }

    public void setReadNegativeStrandFlag(boolean b) {
        mRecord.setReadNegativeStrandFlag(b);
        mNegativeStrandFlag = b;
    }

    /////////////////////////////////////////////////////////////////////
    // *** The following methods are final and cannot be overridden ***//
    /////////////////////////////////////////////////////////////////////

    //public final Object getAttribute(String s) { return mRecord.getAttribute(s); }

    //public final Integer getIntegerAttribute(java.lang.String s) { return mRecord.getIntegerAttribute(s); }

    //public final Short getShortAttribute(String s) { return mRecord.getShortAttribute(s); }

    //public final Byte getByteAttribute(java.lang.String s) { return mRecord.getByteAttribute(s); }

    //public final void setAttribute(String s, Object o) { mRecord.setAttribute(s, o); }

    //public final List<net.sf.samtools.SAMRecord.SAMTagAndValue> getAttributes() { return mRecord.getAttributes(); }


    /////////////////////////////////////////////////////////////////////////////////
    // *** The following methods just call the appropriate method in the record ***//
    /////////////////////////////////////////////////////////////////////////////////

    public String getReadName() { return mRecord.getReadName(); }

    public int getReadNameLength() { return mRecord.getReadNameLength(); }

    public void setReadName(String s) { mRecord.setReadName(s); }

    public byte[] getReadBases() { return mRecord.getReadBases(); }

    public void setReadBases(byte[] bytes) { mRecord.setReadBases(bytes); }

    public int getReadLength() { return mRecord.getReadLength(); }

    public String getBaseQualityString() { return mRecord.getBaseQualityString(); }

    public void setBaseQualityString(java.lang.String s) { mRecord.setBaseQualityString(s); }

    public byte[] getBaseQualities() { return mRecord.getBaseQualities(); }

    public void setBaseQualities(byte[] bytes) { mRecord.setBaseQualities(bytes); }

    public byte[] getOriginalBaseQualities() { return mRecord.getOriginalBaseQualities(); }

    public void setOriginalBaseQualities(byte[] bytes) { mRecord.setOriginalBaseQualities(bytes); }

    public String getReferenceName() { return mRecord.getReferenceName(); }

    public void setReferenceName(String s) { mRecord.setReferenceName(s); }

    public Integer getReferenceIndex() { return mRecord.getReferenceIndex(); }

    public void setReferenceIndex(int i) { mRecord.setReferenceIndex(i); }

    public String getMateReferenceName() { return mRecord.getMateReferenceName(); }

    public void setMateReferenceName(String s) { mRecord.setMateReferenceName(s); }

    public Integer getMateReferenceIndex() { return mRecord.getMateReferenceIndex(); }

    public void setMateReferenceIndex(int i) { mRecord.setMateReferenceIndex(i); }

    public int getAlignmentStart() { return mRecord.getAlignmentStart(); }

    public void setAlignmentStart(int i) { mRecord.setAlignmentStart(i); }

    public int getAlignmentEnd() { return mRecord.getAlignmentEnd(); }

    public int getUnclippedStart() { return mRecord.getUnclippedStart(); }

    public int getUnclippedEnd() { return mRecord.getUnclippedEnd(); }

    public void setAlignmentEnd(int i) { mRecord.setAlignmentEnd(i); }

    public int getMateAlignmentStart() { return mRecord.getMateAlignmentStart(); }

    public void setMateAlignmentStart(int i) { mRecord.setMateAlignmentStart(i); }

    public int getInferredInsertSize() { return mRecord.getInferredInsertSize(); }

    public void setInferredInsertSize(int i) { mRecord.setInferredInsertSize(i); }

    public int getMappingQuality() { return mRecord.getMappingQuality(); }

    public void setMappingQuality(int i) { mRecord.setMappingQuality(i); }

    public String getCigarString() { return mRecord.getCigarString(); }

    public void setCigarString(String s) { mRecord.setCigarString(s); }

    public Cigar getCigar() { return mRecord.getCigar(); }

    public int getCigarLength() { return mRecord.getCigarLength(); }

    public void setCigar(Cigar cigar) { mRecord.setCigar(cigar); }

    public int getFlags() { return mRecord.getFlags(); }

    public void setFlags(int i) { mRecord.setFlags(i); }

    public boolean getReadPairedFlag() { return mRecord.getReadPairedFlag(); }

    public boolean getProperPairFlag() { return mRecord.getProperPairFlag(); }

    public boolean getMateUnmappedFlag() { return mRecord.getMateUnmappedFlag(); }

    public boolean getMateNegativeStrandFlag() { return mRecord.getMateNegativeStrandFlag(); }

    public boolean getFirstOfPairFlag() { return mRecord.getFirstOfPairFlag(); }

    public boolean getSecondOfPairFlag() { return mRecord.getSecondOfPairFlag(); }

    public boolean getNotPrimaryAlignmentFlag() { return mRecord.getNotPrimaryAlignmentFlag(); }

    public boolean getReadFailsVendorQualityCheckFlag() { return mRecord.getReadFailsVendorQualityCheckFlag(); }

    public boolean getDuplicateReadFlag() { return mRecord.getDuplicateReadFlag(); }

    public void setReadPairedFlag(boolean b) { mRecord.setReadPairedFlag(b); }

    public void setProperPairFlag(boolean b) { mRecord.setProperPairFlag(b); }

    public void setMateUnmappedFlag(boolean b) { mRecord.setMateUnmappedFlag(b); }

    public void setMateNegativeStrandFlag(boolean b) { mRecord.setMateNegativeStrandFlag(b); }

    public void setFirstOfPairFlag(boolean b) { mRecord.setFirstOfPairFlag(b); }

    public void setSecondOfPairFlag(boolean b) { mRecord.setSecondOfPairFlag(b); }

    public void setNotPrimaryAlignmentFlag(boolean b) { mRecord.setNotPrimaryAlignmentFlag(b); }

    public void setReadFailsVendorQualityCheckFlag(boolean b) { mRecord.setReadFailsVendorQualityCheckFlag(b); }

    public void setDuplicateReadFlag(boolean b) { mRecord.setDuplicateReadFlag(b); }

    public net.sf.samtools.SAMFileReader.ValidationStringency getValidationStringency() { return mRecord.getValidationStringency(); }

    public void setValidationStringency(net.sf.samtools.SAMFileReader.ValidationStringency validationStringency) { mRecord.setValidationStringency(validationStringency); }

    public SAMFileHeader getHeader() { return mRecord.getHeader(); }

    public void setHeader(SAMFileHeader samFileHeader) { mRecord.setHeader(samFileHeader); }

    public byte[] getVariableBinaryRepresentation() { return mRecord.getVariableBinaryRepresentation(); }

    public int getAttributesBinarySize() { return mRecord.getAttributesBinarySize(); }

    public String format() { return mRecord.format(); }

    public List<AlignmentBlock> getAlignmentBlocks() { return mRecord.getAlignmentBlocks(); }

    public List<SAMValidationError> validateCigar(long l) { return mRecord.validateCigar(l); }

    public boolean equals(Object o) { return mRecord.equals(o); }

    public int hashCode() { return mRecord.hashCode(); }

    public List<SAMValidationError> isValid() { return mRecord.isValid(); }

    public Object clone() throws CloneNotSupportedException { return mRecord.clone(); }

    public String toString() { return mRecord.toString(); }
}
