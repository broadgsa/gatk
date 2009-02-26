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


import edu.mit.broad.sam.util.StringUtil;

import java.util.*;

/**
 * Java binding for a SAM file record.
 */
public class SAMRecord
{
    public static final int UNKNOWN_MAPPING_QUALITY = 255;
    public static final int NO_MAPPING_QUALITY = 0;
    public static final String NO_ALIGNMENT_REFERENCE_NAME = "*";
    public static final String NO_ALIGNMENT_CIGAR = "*";
    public static final int NO_ALIGNMENT_START = 0;
    public static final byte[] NULL_SEQUENCE = "*".getBytes();
    public static final byte[] NULL_QUALS = "*".getBytes();
    private static final int READ_PAIRED_FLAG = 0x1;
    private static final int PROPER_PAIR_FLAG = 0x2;
    private static final int READ_UNMAPPED_FLAG = 0x4;
    private static final int MATE_UNMAPPED_FLAG = 0x8;
    private static final int READ_STRAND_FLAG = 0x10;
    private static final int MATE_STRAND_FLAG = 0x20;
    private static final int FIRST_OF_PAIR_FLAG = 0x40;
    private static final int SECOND_OF_PAIR_FLAG = 0x80;
    private static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    private static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    private static final int DUPLICATE_READ_FLAG = 0x400;


    private String mReadName = null;
    private byte[] mReadBases = NULL_SEQUENCE;
    private byte[] mBaseQualities = NULL_QUALS;
    private String mReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
    private int mAlignmentStart = NO_ALIGNMENT_START;
    private int mMappingQuality = NO_MAPPING_QUALITY;
    private String mCigarString = NO_ALIGNMENT_CIGAR;
    private Cigar mCigar = null;
    private List<AlignmentBlock> mAlignmentBlocks = null;
    private int mFlags = 0;
    private String mMateReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
    private int mMateAlignmentStart = 0;
    private int mInferredInsertSize = 0;
    private Map<String, Object> mAttributes = null;
    private Integer mReferenceIndex = null;
    private Integer mMateReferenceIndex = null;
    private Integer mIndexingBin = null;

    // Optional, but handy for looking of reference indices
    private SAMFileHeader mHeader = null;


    public SAMRecord() {
    }

    public String getReadName() {
        return mReadName;
    }

    /**
     * This method is preferred over getReadName().length(), because for BAMRecord
     * it may be faster.
     * @return length not including a null terminator
     */
    public int getReadNameLength() {
        return mReadName.length();
    }

    public void setReadName(final String value) {
        mReadName = value;
    }

    public String getReadString() {
        return StringUtil.bytesToString(getReadBases());
    }

    public void setReadString(final String value) {
        mReadBases = StringUtil.stringToBytes(value);
    }

    // Read bases, as bytes
    public byte[] getReadBases() {
        return mReadBases;
    }

    public void setReadBases(final byte[] value) {
        mReadBases = value;
    }

    /**
     * This method is preferred over getReadBases().length, because for BAMRecord it may be faster.
     * @return number of bases in the read
     */
    public int getReadLength() {
        return getReadBases().length;
    }

    // Base qualities, encoded as a FASTQ string
    public String getBaseQualityString() {
        return SAMUtils.phredToFastq(getBaseQualities());
    }

    public void setBaseQualityString(final String value) {
        setBaseQualities(SAMUtils.fastqToPhred(value));
    }

    public byte[] getBaseQualities() {
        return mBaseQualities;
    }

    public void setBaseQualities(final byte[] value) {
        mBaseQualities = value;
    }

    public String getReferenceName() {
        return mReferenceName;
    }

    public void setReferenceName(final String value) {
        mReferenceName = value;
        mReferenceIndex = null;
    }

    public Integer getReferenceIndex(final SAMFileHeader header) {
        if (mReferenceIndex == null) {
            if (mReferenceName == null) {
                mReferenceIndex = -1;
            } else if (NO_ALIGNMENT_REFERENCE_NAME.equals(mReferenceName)) {
                mReferenceIndex = -1;
            } else {
                mReferenceIndex = header.getSequenceIndex(mReferenceName);
            }
        }
        return mReferenceIndex;
    }

    public Integer getReferenceIndex() {
        return getReferenceIndex(mHeader);
    }


    public void setReferenceIndex(final int referenceIndex, final SAMFileHeader header) {
        mReferenceIndex = referenceIndex;
        if (mReferenceIndex == -1) {
            mReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
        } else {
            mReferenceName = header.getSequence(referenceIndex).getSequenceName();
        }
    }

    public void setReferenceIndex(final int referenceIndex) {
        setReferenceIndex(referenceIndex, mHeader);
    }


    public String getMateReferenceName() {
        return mMateReferenceName;
    }

    public void setMateReferenceName(final String mateReferenceName) {
        this.mMateReferenceName = mateReferenceName;
        mMateReferenceIndex = null;
    }

    public Integer getMateReferenceIndex(final SAMFileHeader header) {
        if (mMateReferenceIndex == null) {
            if (mMateReferenceName == null) {
                mMateReferenceIndex = -1;
            } else if (NO_ALIGNMENT_REFERENCE_NAME.equals(mMateReferenceName)){
                mMateReferenceIndex = -1;
            } else {
                mMateReferenceIndex = header.getSequenceIndex(mMateReferenceName);
            }
        }
        return mMateReferenceIndex;
    }

    public Integer getMateReferenceIndex() {
        return getMateReferenceIndex(mHeader);
    }

    public void setMateReferenceIndex(final int referenceIndex, final SAMFileHeader header) {
        mMateReferenceIndex = referenceIndex;
        if (mMateReferenceIndex == -1) {
            mMateReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
        } else {
            mMateReferenceName = header.getSequence(referenceIndex).getSequenceName();
        }
    }

    public void setMateReferenceIndex(final int referenceIndex) {
        setMateReferenceIndex(referenceIndex, mHeader);
    }


    public int getAlignmentStart() {
        return mAlignmentStart;
    }

    public void setAlignmentStart(final int value) {
        mAlignmentStart = value;
    }

    public int getAlignmentEnd() {
        final byte[] readBases = getReadBases();
        if (mAlignmentStart == NO_ALIGNMENT_START || Arrays.equals(NULL_SEQUENCE, readBases) || readBases == null) {
            return -1;
        }
        return mAlignmentStart + getCigar().getReferenceLength() - 1;
    }

    /**
     * Returns the alignment start adjusted for clipped bases.  For example if the read
     * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
     * then this method will return 96.
     */
    public int getUnclippedStart() {
        int pos = getAlignmentStart();

        for (final CigarElement cig : getCigar().getCigarElements()) {
            final CigarOperator op = cig.getOperator();
            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                pos -= cig.getLength();
            }
            else {
                break;
            }
        }

        return pos;
    }

    /**
     * Returns the alignment end adjusted for clipped bases.  For example if the read
     * has an alignment end of 100 but the last 7 bases were clipped (hard or soft clipped)
     * then this method will return 107.
     */
    public int getUnclippedEnd() {
        int pos = getAlignmentEnd();
        List<CigarElement> cigs = getCigar().getCigarElements();
        for (int i=cigs.size() - 1; i>=0; --i) {
            final CigarElement cig = cigs.get(i);
            final CigarOperator op = cig.getOperator();

            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                pos += cig.getLength();
            }
            else {
                break;
            }
        }

        return pos;               
    }

    public void setAlignmentEnd(final int value) {
        throw new UnsupportedOperationException("Not supported: setAlignmentEnd");
    }

    public int getMateAlignmentStart() {
        return mMateAlignmentStart;
    }

    public void setMateAlignmentStart(final int mateAlignmentStart) {
        this.mMateAlignmentStart = mateAlignmentStart;
    }

    public int getInferredInsertSize() {
        return mInferredInsertSize;
    }

    public void setInferredInsertSize(final int inferredInsertSize) {
        this.mInferredInsertSize = inferredInsertSize;
    }

    public int getMappingQuality() {
        return mMappingQuality;
    }

    public void setMappingQuality(final int value) {
        mMappingQuality = value;
    }

    public String getCigarString() {
        if (mCigarString == null && getCigar() != null) {
            mCigarString = TextCigarCodec.getSingleton().encode(getCigar());
        }
        return mCigarString;
    }

    public void setCigarString(final String value) {
        mCigarString = value;
        mCigar = null;
    }

    public Cigar getCigar() {
        if (mCigar == null && mCigarString != null) {
            mCigar = TextCigarCodec.getSingleton().decode(mCigarString);
        }
        return mCigar;
    }

    /**
     * This method is preferred over getCigar().getNumElements(), because for BAMRecord it may be faster.
     * @return number of cigar elements (number + operator) in the cigar string
     */
    public int getCigarLength() {
        return getCigar().numCigarElements();
    }

    public void setCigar(final Cigar cigar) {
        this.mCigar = cigar;
        mCigarString = null;
    }

    public int getFlags() {
        return mFlags;
    }

    public void setFlags(final int value) {
        mFlags = value;
    }

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair
     */
    public boolean getReadPairedFlag() {
        return (mFlags & READ_PAIRED_FLAG) != 0;
    }

    private void requireReadPaired() {
        if (!getReadPairedFlag()) {
            throw new IllegalStateException("Inappropriate call if not paired read");
        }
    }

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)
     */
    public boolean getProperPairFlag() {
        requireReadPaired();
        return (mFlags & PROPER_PAIR_FLAG) != 0;
    }

    /**
     * the query sequence itself is unmapped
     */
    public boolean getReadUnmappedFlag() {
        return (mFlags & READ_UNMAPPED_FLAG) != 0;
    }

    /**
     * the mate is unmapped
     */
    public boolean getMateUnmappedFlag() {
        requireReadPaired();
        return (mFlags & MATE_UNMAPPED_FLAG) != 0;
    }

    /**
     * strand of the query (false for forward; true for reverse strand)
     */
    public boolean getReadNegativeStrandFlag() {
        return (mFlags & READ_STRAND_FLAG) != 0;
    }

    /**
     * strand of the mate (false for forward; true for reverse strand)
     */
    public boolean getMateNegativeStrandFlag() {
        requireReadPaired();
        return (mFlags & MATE_STRAND_FLAG) != 0;
    }

    /**
     * the read is the first read in a pair
     */
    public boolean getFirstOfPairFlag() {
        requireReadPaired();
        return (mFlags & FIRST_OF_PAIR_FLAG) != 0;
    }

    /**
     * the read is the second read in a pair
     */
    public boolean getSecondOfPairFlag() {
        requireReadPaired();
        return (mFlags & SECOND_OF_PAIR_FLAG) != 0;
    }

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records)
     */
    public boolean getNotPrimaryAlignmentFlag() {
        return (mFlags & NOT_PRIMARY_ALIGNMENT_FLAG) != 0;
    }

    /**
     * the read fails platform/vendor quality checks
     */
    public boolean getReadFailsVendorQualityCheckFlag() {
        return (mFlags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    /**
     * the read is either a PCR duplicate or an optical duplicate
     */
    public boolean getDuplicateReadFlag() {
        return (mFlags & DUPLICATE_READ_FLAG) != 0;
    }

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair
     */
    public void setReadPairedFlag(final boolean flag) {
        setFlag(flag, READ_PAIRED_FLAG);
    }

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)
     */
    public void setProperPairFlag(final boolean flag) {
        setFlag(flag, PROPER_PAIR_FLAG);
    }

    /**
     * the query sequence itself is unmapped
     */
    public void setReadUmappedFlag(final boolean flag) {
        setFlag(flag, READ_UNMAPPED_FLAG);
    }

    /**
     * the mate is unmapped
     */
    public void setMateUnmappedFlag(final boolean flag) {
        setFlag(flag, MATE_UNMAPPED_FLAG);
    }

    /**
     * strand of the query (false for forward; true for reverse strand)
     */
    public void setReadNegativeStrandFlag(final boolean flag) {
        setFlag(flag, READ_STRAND_FLAG);
    }

    /**
     * strand of the mate (false for forward; true for reverse strand)
     */
    public void setMateNegativeStrandFlag(final boolean flag) {
        setFlag(flag, MATE_STRAND_FLAG);
    }

    /**
     * the read is the first read in a pair
     */
    public void setFirstOfPairFlag(final boolean flag) {
        setFlag(flag, FIRST_OF_PAIR_FLAG);
    }

    /**
     * the read is the second read in a pair
     */
    public void setSecondOfPairFlag(final boolean flag) {
        setFlag(flag, SECOND_OF_PAIR_FLAG);
    }

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records)
     */
    public void setNotPrimaryAlignmentFlag(final boolean flag) {
        setFlag(flag, NOT_PRIMARY_ALIGNMENT_FLAG);
    }

    /**
     * the read fails platform/vendor quality checks
     */
    public void setReadFailsVendorQualityCheckFlag(final boolean flag) {
        setFlag(flag, READ_FAILS_VENDOR_QUALITY_CHECK_FLAG);
    }

    /**
     * the read is either a PCR duplicate or an optical duplicate
     */
    public void setDuplicateReadFlag(final boolean flag) {
        setFlag(flag, DUPLICATE_READ_FLAG);
    }

    private void setFlag(final boolean flag, final int bit) {
        if (flag) {
            mFlags |= bit;
        } else {
            mFlags &= ~bit;
        }
    }

    public Object getAttribute(final String key) {
        if (mAttributes == null) {
            return null;
        }
        return mAttributes.get(key);
    }

    public void setAttribute(final String key, final Object value) {
        if (mAttributes == null) {
            mAttributes = new LinkedHashMap<String, Object>();
        }
        mAttributes.put(key, value);
    }

    public Set<Map.Entry<String, Object>> getAttributes() {
        if (mAttributes == null) {
            return null;
        }
        return mAttributes.entrySet();
    }

    public Integer getIndexingBin() {
        return mIndexingBin;
    }

    public void setIndexingBin(final Integer mIndexingBin) {
        this.mIndexingBin = mIndexingBin;
    }

    public SAMFileHeader getHeader() {
        return mHeader;
    }

    public void setHeader(final SAMFileHeader mHeader) {
        this.mHeader = mHeader;
    }

    /**
     * If this record has a valid binary representation of the variable-length portion of a binary record stored,
     * return that byte array, otherwise return null.  This will never be true for SAMRecords.  It will be true
     * for BAMRecords that have not been eagerDecoded(), and for which none of the data in the variable-length
     * portion has been changed.
     */
    public byte[] getVariableBinaryRepresentation() {
        return null;
    }

    /**
     * Depending on the concrete implementation, the binary file size of attributes may be known without
     * computing them all.
     * @return binary file size of attribute, if known, else -1
     */
    public int getAttributesBinarySize() {
        return -1;
    }

    public String format() {
        final StringBuilder buffer = new StringBuilder();
        addField(buffer, getReadName(), null, null);
        addField(buffer, getFlags(), null, null);
        addField(buffer, getReferenceName(), null, "*");
        addField(buffer, getAlignmentStart(), 0, "*");
        addField(buffer, getMappingQuality(), 0, "0");
        addField(buffer, getCigarString(), null, "*");
        addField(buffer, getMateReferenceName(), null, "*");
        addField(buffer, getMateAlignmentStart(), 0, "*");
        addField(buffer, getInferredInsertSize(), 0, "*");
        addField(buffer, getReadString(), null, "*");
        addField(buffer, getBaseQualityString(), null, "*");
        if (mAttributes != null) {
            for (final Map.Entry<String, Object> entry : getAttributes()) {
                addField(buffer, formatTagValue(entry.getKey(), entry.getValue()));
            }
        }
        return buffer.toString();
    }

    private void addField(final StringBuilder buffer, final Object value, final Object defaultValue, final String defaultString) {
        if (safeEquals(value, defaultValue)) {
            addField(buffer, defaultString);
        } else if (value == null) {
            addField(buffer, "");
        } else {
            addField(buffer, value.toString());
        }
    }

    private void addField(final StringBuilder buffer, final String field) {
        if (buffer.length() > 0) {
            buffer.append('\t');
        }
        buffer.append(field);
    }

    private String formatTagValue(final String key, final Object value) {
        if (value == null || value instanceof String) {
            return key + ":Z:" + value;
        } else if (value instanceof Integer) {
            return key + ":i:" + value;
        } else if (value instanceof Character) {
            return key + ":A:" + value;
        } else if (value instanceof Float) {
            return key + ":f:" + value;
        } else if (value instanceof byte[]) {
            return key + ":H:" + SAMUtils.bytesToHexString((byte[]) value);
        } else {
            throw new RuntimeException("Unexpected value type for key " + key +
                                       ": " + value);
        }
    }

    private boolean safeEquals(final Object o1, final Object o2) {
        if (o1 == o2) {
            return true;
        } else if (o1 == null || o2 == null) {
            return false;
        } else {
            return o1.equals(o2);
        }
    }

    /**
     * Force all lazily-initialized data members to be initialized.  If a subclass overrides this method,
     * typically it should also call  super method.
     */
    protected void eagerDecode() {
        getCigar();
        getCigarString();
    }

    /**
     * Returns blocks of the read sequence that have been aligned directly to the
     * reference sequence. Note that clipped portions of the read and inserted and
     * deleted bases (vs. the reference) are not represented in the alignment blocks.
     */
    public List<AlignmentBlock> getAlignmentBlocks() {
        if (this.mAlignmentBlocks != null) return this.mAlignmentBlocks;

        final Cigar cigar = getCigar();
        if (cigar == null) return Collections.emptyList();


        this.mAlignmentBlocks = new ArrayList<AlignmentBlock>();
        int readBase = 1;
        int refBase  = getAlignmentStart();

        for (final CigarElement e : cigar.getCigarElements()) {
            switch (e.getOperator()) {
                case H : break; // ignore hard clips
                case P : break; // ignore pads
                case S : readBase += e.getLength(); break; // soft clip read bases
                case N : refBase += e.getLength(); break;  // reference skip
                case D : refBase += e.getLength(); break;
                case I : readBase += e.getLength(); break;
                case M :
                    final int length = e.getLength();
                    this.mAlignmentBlocks.add(new AlignmentBlock(readBase, refBase, length));
                    readBase += length;
                    refBase  += length;
                    break;
                default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
            }
        }

        return this.mAlignmentBlocks;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof SAMRecord)) return false;

        final SAMRecord samRecord = (SAMRecord) o;
        eagerDecode();
        samRecord.eagerDecode();

        if (mAlignmentStart != samRecord.mAlignmentStart) return false;
        if (mFlags != samRecord.mFlags) return false;
        if (mInferredInsertSize != samRecord.mInferredInsertSize) return false;
        if (mMappingQuality != samRecord.mMappingQuality) return false;
        if (mMateAlignmentStart != samRecord.mMateAlignmentStart) return false;
        if (mAttributes != null ? !mAttributes.equals(samRecord.mAttributes) : samRecord.mAttributes != null)
            return false;
        if (!Arrays.equals(mBaseQualities, samRecord.mBaseQualities)) return false;
        if (mCigar != null ? !mCigar.equals(samRecord.mCigar) : samRecord.mCigar != null)
            return false;
        if (mIndexingBin != null ? !mIndexingBin.equals(samRecord.mIndexingBin) : samRecord.mIndexingBin != null)
            return false;
        if (mMateReferenceIndex != null ? !mMateReferenceIndex.equals(samRecord.mMateReferenceIndex) : samRecord.mMateReferenceIndex != null)
            return false;
        if (mMateReferenceName != null ? !mMateReferenceName.equals(samRecord.mMateReferenceName) : samRecord.mMateReferenceName != null)
            return false;
        if (!Arrays.equals(mReadBases, samRecord.mReadBases)) return false;
        if (mReadName != null ? !mReadName.equals(samRecord.mReadName) : samRecord.mReadName != null) return false;
        if (mReferenceIndex != null ? !mReferenceIndex.equals(samRecord.mReferenceIndex) : samRecord.mReferenceIndex != null)
            return false;
        if (mReferenceName != null ? !mReferenceName.equals(samRecord.mReferenceName) : samRecord.mReferenceName != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        eagerDecode();
        int result = mReadName != null ? mReadName.hashCode() : 0;
        result = 31 * result + (mReadBases != null ? Arrays.hashCode(mReadBases) : 0);
        result = 31 * result + (mBaseQualities != null ? Arrays.hashCode(mBaseQualities) : 0);
        result = 31 * result + (mReferenceName != null ? mReferenceName.hashCode() : 0);
        result = 31 * result + mAlignmentStart;
        result = 31 * result + mMappingQuality;
        result = 31 * result + (mCigarString != null ? mCigarString.hashCode() : 0);
        result = 31 * result + mFlags;
        result = 31 * result + (mMateReferenceName != null ? mMateReferenceName.hashCode() : 0);
        result = 31 * result + mMateAlignmentStart;
        result = 31 * result + mInferredInsertSize;
        result = 31 * result + (mAttributes != null ? mAttributes.hashCode() : 0);
        result = 31 * result + (mReferenceIndex != null ? mReferenceIndex.hashCode() : 0);
        result = 31 * result + (mMateReferenceIndex != null ? mMateReferenceIndex.hashCode() : 0);
        result = 31 * result + (mIndexingBin != null ? mIndexingBin.hashCode() : 0);
        return result;
    }
}

