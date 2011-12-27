/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.sam;

import com.google.java.contract.Ensures;
import net.sf.samtools.*;
import org.broadinstitute.sting.utils.NGSPlatform;

import java.util.HashMap;
import java.util.Map;

/**
 * @author ebanks, depristo
 * GATKSAMRecord
 *
 * this class extends the samtools BAMRecord class (and SAMRecord) and caches important
 * (and oft-accessed) data that's not already cached by the SAMRecord class
 *
 * IMPORTANT NOTE: Because ReadGroups are not set through the SAMRecord,
 *   if they are ever modified externally then one must also invoke the
 *   setReadGroup() method here to ensure that the cache is kept up-to-date.
 *
 */
public class GATKSAMRecord extends BAMRecord {
    public static final String REDUCED_READ_CONSENSUS_TAG = "RR";

    // the SAMRecord data we're caching
    private String mReadString = null;
    private GATKSAMReadGroupRecord mReadGroup = null;
    private byte[] reducedReadCounts = null;

    // because some values can be null, we don't want to duplicate effort
    private boolean retrievedReadGroup = false;
    private boolean retrievedReduceReadCounts = false;

    // These temporary attributes were added here to make life easier for
    // certain algorithms by providing a way to label or attach arbitrary data to
    // individual GATKSAMRecords.
    // These attributes exist in memory only, and are never written to disk.
    private Map<Object, Object> temporaryAttributes;

    /**
     * HACK TO CREATE GATKSAMRECORD WITH ONLY A HEADER FOR TESTING PURPOSES ONLY
     * @param header
     */
    public GATKSAMRecord(final SAMFileHeader header) {
        this(new SAMRecord(header));
    }

    /**
     * HACK TO CREATE GATKSAMRECORD BASED ONLY A SAMRECORD FOR TESTING PURPOSES ONLY
     * @param read
     */
    public GATKSAMRecord(final SAMRecord read) {
        super(read.getHeader(), read.getMateReferenceIndex(),
                read.getAlignmentStart(),
                read.getReadName() != null ? (short)read.getReadNameLength() : 0,
                (short)read.getMappingQuality(),
                0,
                read.getCigarLength(),
                read.getFlags(),
                read.getReadLength(),
                read.getMateReferenceIndex(),
                read.getMateAlignmentStart(),
                read.getInferredInsertSize(),
                null);
        SAMReadGroupRecord samRG = read.getReadGroup();
        clearAttributes();
        if (samRG != null) {
            GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(samRG);
            setReadGroup(rg);
        }
    }

    public GATKSAMRecord(final SAMFileHeader header,
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
    public GATKSAMReadGroupRecord getReadGroup() {
        if ( !retrievedReadGroup ) {
            SAMReadGroupRecord tempReadGroup = super.getReadGroup();
            mReadGroup = (tempReadGroup == null ? null : new GATKSAMReadGroupRecord(tempReadGroup));
            retrievedReadGroup = true;
        }
        return mReadGroup;
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;

        if (!(o instanceof GATKSAMRecord)) return false;

        // note that we do not consider the GATKSAMRecord internal state at all
        return super.equals(o);
    }

    /**
     * Efficient caching accessor that returns the GATK NGSPlatform of this read
     * @return
     */
    public NGSPlatform getNGSPlatform() {
        return getReadGroup().getNGSPlatform();
    }

    public void setReadGroup( final GATKSAMReadGroupRecord readGroup ) {
        mReadGroup = readGroup;
        retrievedReadGroup = true;
        setAttribute("RG", mReadGroup.getId());       // todo -- this should be standardized, but we don't have access to SAMTagUtils!
    }

    ///////////////////////////////////////////////////////////////////////////////
    // *** ReduceReads functions                                              ***//
    ///////////////////////////////////////////////////////////////////////////////

    public byte[] getReducedReadCounts() {
        if ( ! retrievedReduceReadCounts ) {
            reducedReadCounts = getByteArrayAttribute(REDUCED_READ_CONSENSUS_TAG);
            retrievedReduceReadCounts = true;
        }

        return reducedReadCounts;
    }

    public boolean isReducedRead() {
        return getReducedReadCounts() != null;
    }

    public final byte getReducedCount(final int i) {
        byte firstCount = getReducedReadCounts()[0];
        byte offsetCount = getReducedReadCounts()[i];
        return (i==0) ? firstCount : (byte) Math.min(firstCount + offsetCount, Byte.MAX_VALUE);
    }


    ///////////////////////////////////////////////////////////////////////////////
    // *** GATKSAMRecord specific methods                                     ***//
    ///////////////////////////////////////////////////////////////////////////////


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
     * Checks whether if the read has any bases.
     *
     * Empty reads can be dangerous as it may have no cigar strings, no read names and
     * other missing attributes.
     *
     * @return true if the read has no bases
     */
    public boolean isEmpty() {
        return super.getReadBases() == null || super.getReadLength() == 0;
    }

    /**
     * Clears all attributes except ReadGroup of the read.
     */
    public void simplify () {
        GATKSAMReadGroupRecord rg = getReadGroup();
        this.clearAttributes();
        setReadGroup(rg);
    }

    /**
     * Calculates the reference coordinate for the beginning of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedStart() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped start of the read taking soft clips (but not hard clips) into account
     */
    @Ensures({"result >= getUnclippedStart()", "result <= getUnclippedEnd() || ReadUtils.readIsEntirelyInsertion(this)"})
    public int getSoftStart() {
        int start = this.getUnclippedStart();
        for (CigarElement cigarElement : this.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                start += cigarElement.getLength();
            else
                break;
        }
        return start;
    }

    /**
     * Calculates the reference coordinate for the end of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedStart() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped end of the read taking soft clips (but not hard clips) into account
     */
    @Ensures({"result >= getUnclippedStart()", "result <= getUnclippedEnd() || ReadUtils.readIsEntirelyInsertion(this)"})
    public int getSoftEnd() {
        int stop = this.getUnclippedStart();

        if (ReadUtils.readIsEntirelyInsertion(this))
            return stop;

        int shift = 0;
        CigarOperator lastOperator = null;
        for (CigarElement cigarElement : this.getCigar().getCigarElements()) {
            stop += shift;
            lastOperator = cigarElement.getOperator();
            if (cigarElement.getOperator().consumesReferenceBases() || cigarElement.getOperator() == CigarOperator.SOFT_CLIP || cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                shift = cigarElement.getLength();
            else
                shift = 0;
        }
        return (lastOperator == CigarOperator.HARD_CLIP) ? stop-1 : stop+shift-1 ;
    }



}
