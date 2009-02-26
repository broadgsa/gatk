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

import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.StringUtil;

import java.io.ByteArrayInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;


/**
 * Wrapper class for binary BAM records.
 * Delays unpacking all data binary until requested.
 */
class BAMRecord
    extends SAMRecord
{
    private static final int READ_NAME_OFFSET = 0;

    private byte[] mRestOfBinaryData = null;
    private int mReadLength = 0;
    private final short mReadNameLength;
    private final int mCigarLen;
    private boolean mAttributesDecoded = false;
    private boolean mCigarDecoded = false;

    /**
     * If any of the properties set from mRestOfBinaryData have been overridden by calls to setters,
     * this is set to true, indicating that mRestOfBinaryData cannot be used to write this record to disk.
     */
    private boolean mBinaryDataStale;

    BAMRecord(final SAMFileHeader header, final int referenceID, final int coordinate, final short readNameLength, final short mappingQuality,
              final int indexingBin, final int cigarLen, final int flags, final int readLen, final int mateReferenceID, final int mateCoordinate, final int insertSize,
              final byte[] restOfData) {
        setReferenceIndex(referenceID, header);
        setAlignmentStart(coordinate);
        mReadNameLength = readNameLength;
        setMappingQuality(mappingQuality);
        setIndexingBin(indexingBin);
        mCigarLen = cigarLen;
        setFlags(flags);
        mReadLength = readLen;
        setMateReferenceIndex(mateReferenceID, header);
        setMateAlignmentStart(mateCoordinate);
        setInferredInsertSize(insertSize);
        mRestOfBinaryData = restOfData;

        // Set these to null in order to mark them as being candidates for lazy initialization.
        // If this is not done, they will have non-null defaults.
        super.setReadName(null);
        super.setCigarString(null);
        super.setReadBases(null);
        super.setBaseQualities(null);

        // Mark the binary block as being valid for writing back out to disk
        mBinaryDataStale = false;
    }

    protected void eagerDecode() {
        // Force all the lazily-initialized attributes to be decoded.
        getReadName();
        getCigar();
        getReadBases();
        getBaseQualities();
        getAttributes();
        super.eagerDecode();
        mRestOfBinaryData = null;
    }

    /**
     * If this record has a valid binary representation of the variable-length portion of a binary record stored,
     * return that byte array, otherwise return null.  This will never be true for SAMRecords.  It will be true
     * for BAMRecords that have not been eagerDecoded(), and for which none of the data in the variable-length
     * portion has been changed.
     */
    @Override
    public byte[] getVariableBinaryRepresentation() {
        if (mBinaryDataStale) {
            return null;
        }
        // This may have been set to null by eagerDecode()
        return mRestOfBinaryData;
    }

    /**
     * Depending on the concrete implementation, the binary file size of attributes may be known without
     * computing them all.
     *
     * @return binary file size of attribute, if known, else -1
     */
    @Override
    public int getAttributesBinarySize() {
        if (mBinaryDataStale || mRestOfBinaryData == null) {
            return -1;
        }
        final int tagsOffset = readNameSize() + cigarSize() + basesSize() + qualsSize();
        return mRestOfBinaryData.length - tagsOffset;
    }

    @Override
    public void setReadName(final String value) {
        super.setReadName(value);
        mBinaryDataStale = true;
    }

    @Override
    public void setCigar(final Cigar cigar) {
        super.setCigar(cigar);
        mBinaryDataStale = true;
    }

    @Override
    public void setReadBases(final byte[] value) {
        super.setReadBases(value);
        mBinaryDataStale = true;
    }

    @Override
    public void setBaseQualities(final byte[] value) {
        super.setBaseQualities(value);
        mBinaryDataStale = true;
    }

    @Override
    public void setAttribute(final String key, final Object value) {
        // populate all the attributes from the binary block before overwriting one
        getAttributes();
        super.setAttribute(key, value);
        mBinaryDataStale = true;
    }

    /**
     * Avoids decoding binary block to get read length
     */
    @Override
    public int getReadLength() {
        return mReadLength;
    }

    @Override
    public String getReadName() {
        String result = super.getReadName();
        if (mRestOfBinaryData != null && result == null) {
            result = decodeReadName();
            super.setReadName(result);
        }
        return result;
    }

    /**
     * Do not include null terminator
     */
    @Override
    public int getReadNameLength() {
        return mReadNameLength - 1;
    }

    @Override
    public Cigar getCigar() {
        if (mRestOfBinaryData != null && !mCigarDecoded) {
            final int cigarOffset = readNameSize();
            final ByteBuffer byteBuffer  = ByteBuffer.wrap(mRestOfBinaryData, cigarOffset, cigarSize());
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
            super.setCigar(BinaryCigarCodec.getSingleton().decode(byteBuffer));
            mCigarDecoded = true;
        }
        return super.getCigar();
    }

    @Override
    public int getCigarLength() {
        return mCigarLen;
    }

    @Override
    public byte[] getReadBases() {
        byte[] result = super.getReadBases();
        if (mRestOfBinaryData != null && result == null && mReadLength > 0) {
            result = decodeReadBases();
            super.setReadBases(result);
        }
        return result;
    }

    @Override
    public byte[] getBaseQualities() {
        byte[] ret = super.getBaseQualities();
        if (mRestOfBinaryData != null && ret == null && mReadLength > 0) {
            ret = decodeBaseQualities();
            super.setBaseQualities(ret);
        }
        return ret;
    }

    @Override
    public Object getAttribute(final String key) {
        if (!mAttributesDecoded) {
            decodeAttributes();
        }
        return super.getAttribute(key);
    }

    @Override
    public Set<Map.Entry<String, Object>> getAttributes() {
        if (!mAttributesDecoded) {
            decodeAttributes();
        }
        return super.getAttributes();
    }

    private void decodeAttributes() {
        if (mAttributesDecoded) {
            return;
        }
        mAttributesDecoded = true;
        final Map<String,Object> attributes = new LinkedHashMap<String, Object>();
        final int tagsOffset = readNameSize() + cigarSize() + basesSize() + qualsSize();
        final int tagsSize = mRestOfBinaryData.length - tagsOffset;
        final BinaryCodec byteBufferCodec = new BinaryCodec(new ByteArrayInputStream(mRestOfBinaryData, tagsOffset, tagsSize));
        new BinaryTagCodec(byteBufferCodec).readTags(attributes);
        for (final Map.Entry<String, Object> entry : attributes.entrySet()) {
            super.setAttribute(entry.getKey(), entry.getValue());
        }
    }

    private byte[] decodeBaseQualities() {
        if (mReadLength == 0) {
            return null;
        }
        final int qualsOffset = readNameSize() + cigarSize() + basesSize();
        final byte[] ret = new byte[qualsSize()];
        System.arraycopy(mRestOfBinaryData, qualsOffset, ret, 0, qualsSize());
        return ret;
    }

    private String decodeReadName() {
        // Don't include terminating null
        return StringUtil.bytesToString(mRestOfBinaryData, READ_NAME_OFFSET, mReadNameLength-1);
    }

    private byte[] decodeReadBases() {
        if (mReadLength == 0) {
            return null;
        }
        final int basesOffset = readNameSize() + cigarSize();
        return SAMUtils.compressedBasesToBytes(mReadLength, mRestOfBinaryData, basesOffset);
    }

    /* methods for computing size of variably-sizes elements */

    private int readNameSize() {
        return mReadNameLength;
    }

    private int cigarSize() {
        return mCigarLen * 4;
    }

    private int basesSize() {
        return (mReadLength + 1)/2;
    }

    private int qualsSize() {
        return mReadLength;
    }
}
