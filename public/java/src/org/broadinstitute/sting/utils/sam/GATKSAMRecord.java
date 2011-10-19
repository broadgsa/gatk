package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.*;

import java.util.HashMap;
import java.util.Map;

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
public class GATKSamRecord extends BAMRecord {
    // the SAMRecord data we're caching
    private String mReadString = null;
    private SAMReadGroupRecord mReadGroup = null;

    // because some values can be null, we don't want to duplicate effort
    private boolean retrievedReadGroup = false;

    // These temporary attributes were added here to make life easier for
    // certain algorithms by providing a way to label or attach arbitrary data to
    // individual GATKSAMRecords.
    // These attributes exist in memory only, and are never written to disk.
    private Map<Object, Object> temporaryAttributes;

    /**
     * HACK TO CREATE GATKSAMRECORD WITH ONLY A HEADER FOR TESTING PURPOSES ONLY
     * @param header
     */
    public GATKSamRecord(final SAMFileHeader header) {
        super(header, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, SAMRecord.NO_ALIGNMENT_START,
                (short)0, (short)255, 0, 1, 0, 1, 0, 0, 0, null);
    }


    public GATKSamRecord(final SAMFileHeader header,
                         final int referenceSequenceIndex,
                         final int alignmentStart,
                         final short readNameLength,
                         final short mappingQuality,
                         final int indexingBin,
                         final int cigarLen,
                         final int flags,
                         final int readLen,
                         final int mateReferenceSequenceIndex,
                         final int mateAlignmentStart,
                         final int insertSize,
                         final byte[] variableLengthBlock) {
        super(header, referenceSequenceIndex, alignmentStart, readNameLength, mappingQuality, indexingBin, cigarLen,
                flags, readLen, mateReferenceSequenceIndex, mateAlignmentStart, insertSize, variableLengthBlock);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // *** The following methods are overloaded to cache the appropriate data ***//
    ///////////////////////////////////////////////////////////////////////////////

    @Override
    public String getReadString() {
        if ( mReadString == null )
            mReadString = super.getReadString();
        return mReadString;
    }

    @Override
    public void setReadString(String s) {
        super.setReadString(s);
        mReadString = s;
    }

    @Override
    public SAMReadGroupRecord getReadGroup() {
        if ( !retrievedReadGroup ) {
            SAMReadGroupRecord tempReadGroup = super.getReadGroup();
            mReadGroup = (tempReadGroup == null ? tempReadGroup : new GATKSAMReadGroupRecord(tempReadGroup));
            retrievedReadGroup = true;
        }
        return mReadGroup;
    }

    public void setReadGroup(SAMReadGroupRecord record) {
        mReadGroup = record;
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;

        // note -- this forbids a GATKSAMRecord being equal to its underlying SAMRecord
        if (!(o instanceof GATKSamRecord)) return false;

        // note that we do not consider the GATKSAMRecord internal state at all
        return super.equals(o);
    }
}
