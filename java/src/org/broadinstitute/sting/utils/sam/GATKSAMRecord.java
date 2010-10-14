package org.broadinstitute.sting.utils.sam;

import java.lang.reflect.Method;
import java.util.*;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

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
 * 13 Oct 2010 - mhanna - this class is fundamentally flawed: it uses a decorator
 *                        pattern to wrap a heavyweight object, which can lead
 *                        to heinous side effects if the wrapping is not carefully
 *                        done.  Hopefully SAMRecord will become an interface and
 *                        this will eventually be fixed.
 */
public class GATKSAMRecord extends SAMRecord {

    // the underlying SAMRecord which we are wrapping
    private final SAMRecord mRecord;

    // the SAMRecord data we're caching
    private String mReadString = null;
    private SAMReadGroupRecord mReadGroup = null;
    private boolean mNegativeStrandFlag;
    private boolean mUnmappedFlag;
    private Boolean mSecondOfPairFlag = null;

    // because some values can be null, we don't want to duplicate effort
    private boolean retrievedReadGroup = false;

    // These temporary attributes were added here to make life easier for
    // certain algorithms by providing a way to label or attach arbitrary data to
    // individual GATKSAMRecords.
    // These attributes exist in memory only, and are never written to disk.
    private Map<Object, Object> temporaryAttributes;

    // A bitset which represents the bases of the read.  If a bit is set, then
    // the base is good; otherwise it is a bad base (as defined by the setter).
    // TODO: this is a temporary hack.  If it works, clean it up.
    private BitSet mBitSet = null;

    public GATKSAMRecord(SAMRecord record, boolean useOriginalBaseQualities) {
        super(null); // it doesn't matter - this isn't used
        if ( record == null )
            throw new IllegalArgumentException("The SAMRecord argument cannot be null");
        mRecord = record;

        mNegativeStrandFlag = mRecord.getReadNegativeStrandFlag();
        mUnmappedFlag = mRecord.getReadUnmappedFlag();

        // because attribute methods are declared to be final (and we can't overload them),
        // we need to actually set all of the attributes here
        List<SAMTagAndValue> attributes = record.getAttributes();
        for ( SAMTagAndValue attribute : attributes )
            setAttribute(attribute.tag, attribute.value);

        // if we are using original quals, set them now if they are present in the record
        if ( useOriginalBaseQualities ) {
            byte[] originalQuals = mRecord.getOriginalBaseQualities();
            if ( originalQuals != null )
                mRecord.setBaseQualities(originalQuals);
        }

        // sanity check that the lengths of the base and quality strings are equal
        if ( getBaseQualities().length  != getReadLength() )
            throw new UserException.MalformedBam(this, String.format("Error: the number of base qualities does not match the number of bases in %s (and the GATK does not currently support '*' for the quals)", mRecord.getReadName()));
    }

    public void setGoodBases(GATKSAMRecordFilter filter, boolean abortIfAlreadySet) {
        if ( mBitSet == null || !abortIfAlreadySet )
            mBitSet = filter.getGoodBases(this);
    }

    public boolean isGoodBase(int index) {
        return ( mBitSet == null || mBitSet.size() <= index ? true : mBitSet.get(index));
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
        return mUnmappedFlag;
    }

    public void setReadUmappedFlag(boolean b) {
        mRecord.setReadUmappedFlag(b);
        mUnmappedFlag = b;
    }

    public boolean getReadNegativeStrandFlag() {
        return mNegativeStrandFlag;
    }

    public void setReadNegativeStrandFlag(boolean b) {
        mRecord.setReadNegativeStrandFlag(b);
        mNegativeStrandFlag = b;
    }

    public boolean getSecondOfPairFlag() {
        if( mSecondOfPairFlag == null ) {
            //not done in constructor because this method can't be called for
            //all SAMRecords.
            mSecondOfPairFlag = mRecord.getSecondOfPairFlag();
        }
        return mSecondOfPairFlag;
    }

    public void setSecondOfPairFlag(boolean b) {
        mRecord.setSecondOfPairFlag(b);
        mSecondOfPairFlag = b;
    }

    /**
     * Checks whether an attribute has been set for the given key.
     *
     * Temporary attributes provide a way to label or attach arbitrary data to
     * individual GATKSAMRecords. These attributes exist in memory only,
     * and are never written to disk.
     *
     * @param key key
     * @return True if an attribute has been set for this key.
     */
    public boolean containsTemporaryAttribute(Object key) {
        if(temporaryAttributes != null) {
            return temporaryAttributes.containsKey(key);
        }
        return false;
    }

    /**
     * Sets the key to the given value, replacing any previous value. The previous
     * value is returned.
     *
     * Temporary attributes provide a way to label or attach arbitrary data to
     * individual GATKSAMRecords. These attributes exist in memory only,
     * and are never written to disk.
     *
     * @param key    key
     * @param value  value
     * @return attribute
     */
    public Object setTemporaryAttribute(Object key, Object value) {
        if(temporaryAttributes == null) {
            temporaryAttributes = new HashMap<Object, Object>();
        }
        return temporaryAttributes.put(key, value);
    }

    /**
     * Looks up the value associated with the given key.
     *
     * Temporary attributes provide a way to label or attach arbitrary data to
     * individual GATKSAMRecords. These attributes exist in memory only,
     * and are never written to disk.
     *
     * @param key key
     * @return The value, or null.
     */
    public Object getTemporaryAttribute(Object key) {
        if(temporaryAttributes != null) {
            return temporaryAttributes.get(key);
        }
        return null;
    }

    /**
     * Removes the attribute that has the given key.
     *
     * Temporary attributes provide a way to label or attach arbitrary data to
     * individual GATKSAMRecords. These attributes exist in memory only,
     * and are never written to disk.
     *
     * @param key key
     * @return The value that was associated with this key, or null.
     */
    public Object removeTemporaryAttribute(Object key) {
         if(temporaryAttributes != null) {
             return temporaryAttributes.remove(key);
         }
         return null;
    }

    /////////////////////////////////////////////////////////////////////////////////
    // *** The following methods just call the appropriate method in the record ***//
    /////////////////////////////////////////////////////////////////////////////////

    public String getReadName() { return mRecord.getReadName(); }

    public int getReadNameLength() { return mRecord.getReadNameLength(); }

    public void setReadName(String s) { mRecord.setReadName(s); }

    public byte[] getReadBases() { return mRecord.getReadBases(); }

    public void setReadBases(byte[] bytes) { mRecord.setReadBases(bytes); }

    public int getReadLength() { return mRecord.getReadLength(); }

    public byte[] getBaseQualities() { return mRecord.getBaseQualities(); }

    public void setBaseQualities(byte[] bytes) { mRecord.setBaseQualities(bytes); }

    public String getBaseQualityString() { return mRecord.getBaseQualityString(); }

    public void setBaseQualityString(String s) { mRecord.setBaseQualityString(s); }

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

    public boolean getNotPrimaryAlignmentFlag() { return mRecord.getNotPrimaryAlignmentFlag(); }

    public boolean getReadFailsVendorQualityCheckFlag() { return mRecord.getReadFailsVendorQualityCheckFlag(); }

    public boolean getDuplicateReadFlag() { return mRecord.getDuplicateReadFlag(); }

    public void setReadPairedFlag(boolean b) { mRecord.setReadPairedFlag(b); }

    public void setProperPairFlag(boolean b) { mRecord.setProperPairFlag(b); }

    public void setMateUnmappedFlag(boolean b) { mRecord.setMateUnmappedFlag(b); }

    public void setMateNegativeStrandFlag(boolean b) { mRecord.setMateNegativeStrandFlag(b); }

    public void setFirstOfPairFlag(boolean b) { mRecord.setFirstOfPairFlag(b); }

    public void setNotPrimaryAlignmentFlag(boolean b) { mRecord.setNotPrimaryAlignmentFlag(b); }

    public void setReadFailsVendorQualityCheckFlag(boolean b) { mRecord.setReadFailsVendorQualityCheckFlag(b); }

    public void setDuplicateReadFlag(boolean b) { mRecord.setDuplicateReadFlag(b); }

    public net.sf.samtools.SAMFileReader.ValidationStringency getValidationStringency() { return mRecord.getValidationStringency(); }

    public void setValidationStringency(net.sf.samtools.SAMFileReader.ValidationStringency validationStringency) { mRecord.setValidationStringency(validationStringency); }

    public Object getAttribute(final String tag) { return mRecord.getAttribute(tag); }

    public Integer getIntegerAttribute(final String tag) { return mRecord.getIntegerAttribute(tag); }

    public Short getShortAttribute(final String tag) { return mRecord.getShortAttribute(tag); }

    public Byte getByteAttribute(final String tag) { return mRecord.getByteAttribute(tag); }

    public String getStringAttribute(final String tag) { return mRecord.getStringAttribute(tag); }

    public Character getCharacterAttribute(final String tag) { return mRecord.getCharacterAttribute(tag); }

    public Float getFloatAttribute(final String tag) { return mRecord.getFloatAttribute(tag); }

    public byte[] getByteArrayAttribute(final String tag) { return mRecord.getByteArrayAttribute(tag); }

    protected Object getAttribute(final short tag) {
        Object attribute;
        try {
            Method method = mRecord.getClass().getDeclaredMethod("getAttribute",Short.TYPE);
            method.setAccessible(true);
            attribute = method.invoke(mRecord,tag);
        }
        catch(Exception ex) {
            throw new ReviewedStingException("Unable to invoke getAttribute method",ex);
        }
        return attribute;
    }

    public void setAttribute(final String tag, final Object value) { mRecord.setAttribute(tag,value); }

    protected void setAttribute(final short tag, final Object value) {
        try {
            Method method = mRecord.getClass().getDeclaredMethod("setAttribute",Short.TYPE,Object.class);
            method.setAccessible(true);
            method.invoke(mRecord,tag,value);
        }
        catch(Exception ex) {
            throw new ReviewedStingException("Unable to invoke setAttribute method",ex);
        }
    }

    public void clearAttributes() { mRecord.clearAttributes(); }

    protected void setAttributes(final SAMBinaryTagAndValue attributes) {
        try {
            Method method = mRecord.getClass().getDeclaredMethod("setAttributes",SAMBinaryTagAndValue.class);
            method.setAccessible(true);
            method.invoke(mRecord,attributes);
        }
        catch(Exception ex) {
            throw new ReviewedStingException("Unable to invoke setAttributes method",ex);
        }
   }

    protected SAMBinaryTagAndValue getBinaryAttributes() {
        SAMBinaryTagAndValue binaryAttributes;
        try {
            Method method = mRecord.getClass().getDeclaredMethod("getBinaryAttributes");
            method.setAccessible(true);
            binaryAttributes = (SAMBinaryTagAndValue)method.invoke(mRecord);
        }
        catch(Exception ex) {
            throw new ReviewedStingException("Unable to invoke getBinaryAttributes method",ex);
        }
        return binaryAttributes;
    }

    public List<SAMTagAndValue> getAttributes() { return mRecord.getAttributes(); }

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

    public SAMFileSource getFileSource() { return mRecord.getFileSource(); }
}
