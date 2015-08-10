/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.sam;

import com.google.java.contract.Ensures;
import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.NGSPlatform;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.recalibration.EventType;

import java.util.*;

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
 * WARNING -- GATKSAMRecords cache several values (that are expensive to compute)
 * that depending on the inferred insert size and alignment starts and stops of this read and its mate.
 * Changing these values in any way will invalidate the cached value. However, we do not monitor those setter
 * functions, so modifying a GATKSAMRecord in any way may result in stale cached values.
 */
public class GATKSAMRecord extends SAMRecord implements Cloneable {
    // Base Quality Score Recalibrator specific attribute tags
    public static final String BQSR_BASE_INSERTION_QUALITIES = "BI";                // base qualities for insertions
    public static final String BQSR_BASE_DELETION_QUALITIES = "BD";                 // base qualities for deletions

    /**
     * The default quality score for an insertion or deletion, if
     * none are provided for this read.
     */
    public static final byte DEFAULT_INSERTION_DELETION_QUAL = (byte)45;

    // the SAMRecord data we're caching
    private String mReadString = null;
    private GATKSAMReadGroupRecord mReadGroup = null;
    private final static int UNINITIALIZED = -1;
    private int softStart = UNINITIALIZED;
    private int softEnd = UNINITIALIZED;
    private Integer adapterBoundary = null;

    private boolean isStrandlessRead = false;

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
    public GATKSAMRecord(final SAMFileHeader header) {
        this(new SAMRecord(header));
    }

    /**
     * HACK TO CREATE GATKSAMRECORD BASED ONLY A SAMRECORD FOR TESTING PURPOSES ONLY
     * @param read
     */
    public GATKSAMRecord(final SAMRecord read) {
        super(read.getHeader());
        super.setReferenceIndex(read.getReferenceIndex());
        super.setAlignmentStart(read.getAlignmentStart());
        super.setReadName(read.getReadName());
        super.setMappingQuality(read.getMappingQuality());
        // indexing bin done below
        super.setCigar(read.getCigar());
        super.setFlags(read.getFlags());
        super.setMateReferenceIndex(read.getMateReferenceIndex());
        super.setMateAlignmentStart(read.getMateAlignmentStart());
        super.setInferredInsertSize(read.getInferredInsertSize());
        SAMReadGroupRecord samRG = read.getReadGroup();
        SAMBinaryTagAndValue samAttr = GATKBin.getReadBinaryAttributes(read);
        if (samAttr == null) {
            clearAttributes();
        } else {
            setAttributes(samAttr);
        }
        if (samRG != null) {
            GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(samRG);
            setReadGroup(rg);
        }

        super.setFileSource(read.getFileSource());
        super.setReadName(read.getReadName());
        super.setCigarString(read.getCigarString());
        super.setReadBases(read.getReadBases());
        super.setBaseQualities(read.getBaseQualities());
        // From SAMRecord constructor: Do this after the above because setCigarString will clear it.
        GATKBin.setReadIndexingBin(this, GATKBin.getReadIndexingBin(read));
    }

    public static GATKSAMRecord createRandomRead(int length) {
        List<CigarElement> cigarElements = new LinkedList<>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        return ArtificialSAMUtils.createArtificialRead(cigar);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // *** support for reads without meaningful strand information            ***//
    ///////////////////////////////////////////////////////////////////////////////

    /**
     * Does this read have a meaningful strandedness value?
     *
     * Some advanced types of reads, such as reads coming from merged fragments,
     * don't have meaningful strandedness values, as they are composites of multiple
     * other reads.  Strandless reads need to be handled specially by code that cares about
     * stranded information, such as FS.
     *
     * @return true if this read doesn't have meaningful strand information
     */
    public boolean isStrandless() {
        return isStrandlessRead;
    }

    /**
     * Set the strandless state of this read to isStrandless
     * @param isStrandless true if this read doesn't have a meaningful strandedness value
     */
    public void setIsStrandless(final boolean isStrandless) {
        this.isStrandlessRead = isStrandless;
    }

    @Override
    public boolean getReadNegativeStrandFlag() {
        return ! isStrandless() && super.getReadNegativeStrandFlag();
    }

    @Override
    public void setReadNegativeStrandFlag(final boolean flag) {
        if ( isStrandless() )
            throw new IllegalStateException("Cannot set the strand of a strandless read");
        super.setReadNegativeStrandFlag(flag);
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

    /**
     * Get the GATKSAMReadGroupRecord of this read
     * @return a non-null GATKSAMReadGroupRecord
     */
    @Override
    public GATKSAMReadGroupRecord getReadGroup() {
        if ( ! retrievedReadGroup ) {
            final SAMReadGroupRecord rg = super.getReadGroup();

            // three cases: rg may be null (no rg, rg may already be a GATKSAMReadGroupRecord, or it may be
            // a regular SAMReadGroupRecord in which case we have to make it a GATKSAMReadGroupRecord
            if ( rg == null )
                mReadGroup = null;
            else if ( rg instanceof GATKSAMReadGroupRecord )
                mReadGroup = (GATKSAMReadGroupRecord)rg;
            else
                mReadGroup = new GATKSAMReadGroupRecord(rg);

            retrievedReadGroup = true;
        }
        return mReadGroup;
    }

    public void setReadGroup( final GATKSAMReadGroupRecord readGroup ) {
        mReadGroup = readGroup;
        retrievedReadGroup = true;
        setAttribute("RG", mReadGroup.getId()); // todo -- this should be standardized, but we don't have access to SAMTagUtils!
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
     * Setters and Accessors for base insertion and base deletion quality scores
     */
    public void setBaseQualities( final byte[] quals, final EventType errorModel ) {
        switch( errorModel ) {
            case BASE_SUBSTITUTION:
                setBaseQualities(quals);
                break;
            case BASE_INSERTION:
                setAttribute( GATKSAMRecord.BQSR_BASE_INSERTION_QUALITIES, quals == null ? null : SAMUtils.phredToFastq(quals) );
                break;
            case BASE_DELETION:
                setAttribute( GATKSAMRecord.BQSR_BASE_DELETION_QUALITIES, quals == null ? null : SAMUtils.phredToFastq(quals) );
                break;
            default:
                throw new ReviewedGATKException("Unrecognized Base Recalibration type: " + errorModel );
        }
    }

    public byte[] getBaseQualities( final EventType errorModel ) {
        switch( errorModel ) {
            case BASE_SUBSTITUTION:
                return getBaseQualities();
            case BASE_INSERTION:
                return getBaseInsertionQualities();
            case BASE_DELETION:
                return getBaseDeletionQualities();
            default:
                throw new ReviewedGATKException("Unrecognized Base Recalibration type: " + errorModel );
        }
    }

    /**
     * @return whether or not this read has base insertion or deletion qualities (one of the two is sufficient to return true)
     */
    public boolean hasBaseIndelQualities() {
        return getAttribute( BQSR_BASE_INSERTION_QUALITIES ) != null || getAttribute( BQSR_BASE_DELETION_QUALITIES ) != null;
    }

    /**
     * @return the base deletion quality or null if read doesn't have one
     */
    public byte[] getExistingBaseInsertionQualities() {
        return SAMUtils.fastqToPhred( getStringAttribute(BQSR_BASE_INSERTION_QUALITIES));
    }

    /**
     * @return the base deletion quality or null if read doesn't have one
     */
    public byte[] getExistingBaseDeletionQualities() {
        return SAMUtils.fastqToPhred( getStringAttribute(BQSR_BASE_DELETION_QUALITIES));
    }

    /**
     * Default utility to query the base insertion quality of a read. If the read doesn't have one, it creates an array of default qualities (currently Q45)
     * and assigns it to the read.
     *
     * @return the base insertion quality array
     */
    public byte[] getBaseInsertionQualities() {
        byte [] quals = getExistingBaseInsertionQualities();
        if( quals == null ) {
            quals = new byte[getBaseQualities().length];
            Arrays.fill(quals, DEFAULT_INSERTION_DELETION_QUAL); // Some day in the future when base insertion and base deletion quals exist the samtools API will
                                           // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
        }
        return quals;
    }

    /**
     * Default utility to query the base deletion quality of a read. If the read doesn't have one, it creates an array of default qualities (currently Q45)
     * and assigns it to the read.
     *
     * @return the base deletion quality array
     */
    public byte[] getBaseDeletionQualities() {
        byte[] quals = getExistingBaseDeletionQualities();
        if( quals == null ) {
            quals = new byte[getBaseQualities().length];
            Arrays.fill(quals, DEFAULT_INSERTION_DELETION_QUAL);  // Some day in the future when base insertion and base deletion quals exist the samtools API will
                                            // be updated and the original quals will be pulled here, but for now we assume the original quality is a flat Q45
        }
        return quals;
    }

    /**
     * Efficient caching accessor that returns the GATK NGSPlatform of this read
     * @return
     */
    public NGSPlatform getNGSPlatform() {
        return getReadGroup().getNGSPlatform();
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
        return temporaryAttributes != null && temporaryAttributes.containsKey(key);
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
            temporaryAttributes = new HashMap<>();
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
    public GATKSAMRecord simplify () {
        GATKSAMReadGroupRecord rg = getReadGroup(); // save the read group information
        byte[] insQuals = (this.getAttribute(BQSR_BASE_INSERTION_QUALITIES) == null) ? null : getBaseInsertionQualities();
        byte[] delQuals = (this.getAttribute(BQSR_BASE_DELETION_QUALITIES)  == null) ? null : getBaseDeletionQualities();
        this.clearAttributes(); // clear all attributes from the read
        this.setReadGroup(rg); // restore read group
        if (insQuals != null)
           this.setBaseQualities(insQuals, EventType.BASE_INSERTION); // restore base insertion if we had any
        if (delQuals != null)
            this.setBaseQualities(delQuals, EventType.BASE_DELETION); // restore base deletion if we had any
        return this;
    }

    /**
     * Calculates the reference coordinate for the beginning of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedStart() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped start of the read taking soft clips (but not hard clips) into account
     */
    public int getSoftStart() {
        if ( softStart == UNINITIALIZED ) {
            softStart = getAlignmentStart();
            for (final CigarElement cig : getCigar().getCigarElements()) {
                final CigarOperator op = cig.getOperator();

                if (op == CigarOperator.SOFT_CLIP)
                    softStart -= cig.getLength();
                else if (op != CigarOperator.HARD_CLIP)
                    break;
            }
        }
        return softStart;
    }

    /**
     * Calculates the reference coordinate for the end of the read taking into account soft clips but not hard clips.
     *
     * Note: getUnclippedEnd() adds soft and hard clips, this function only adds soft clips.
     *
     * @return the unclipped end of the read taking soft clips (but not hard clips) into account
     */
    public int getSoftEnd() {
        if ( softEnd == UNINITIALIZED ) {
            boolean foundAlignedBase = false;
            softEnd = getAlignmentEnd();
            final List<CigarElement> cigs = getCigar().getCigarElements();
            for (int i = cigs.size() - 1; i >= 0; --i) {
                final CigarElement cig = cigs.get(i);
                final CigarOperator op = cig.getOperator();

                if (op == CigarOperator.SOFT_CLIP) // assumes the soft clip that we found is at the end of the aligned read
                    softEnd += cig.getLength();
                else if (op != CigarOperator.HARD_CLIP) {
                    foundAlignedBase = true;
                    break;
                }
            }
            if( !foundAlignedBase ) { // for example 64H14S, the soft end is actually the same as the alignment end
                softEnd = getAlignmentEnd();
            }
        }

        return softEnd;
    }

    /**
     * If the read is hard clipped, the soft start and end will change. You can set manually or just reset the cache
     * so that the next call to getSoftStart/End will recalculate it lazily.
     */
    public void resetSoftStartAndEnd() {
        softStart = -1;
        softEnd = -1;
    }

    /**
     * If the read is hard clipped, the soft start and end will change. You can set manually or just reset the cache
     * so that the next call to getSoftStart/End will recalculate it lazily.
     */
    public void resetSoftStartAndEnd(int softStart, int softEnd) {
        this.softStart = softStart;
        this.softEnd = softEnd;
    }

    /**
     * Determines the original alignment start of a previously clipped read.
     * 
     * This is useful for reads that have been trimmed to a variant region and lost the information of it's original alignment end
     * 
     * @return the alignment start of a read before it was clipped
     */
    public int getOriginalAlignmentStart() {
        return getUnclippedStart();
    }

    /**
     * Determines the original alignment end of a previously clipped read.
     *
     * This is useful for reads that have been trimmed to a variant region and lost the information of it's original alignment end
     * 
     * @return the alignment end of a read before it was clipped
     */
    public int getOriginalAlignmentEnd() {
        return getUnclippedEnd();
    }

    /**
     * Creates an empty GATKSAMRecord with the read's header, read group and mate
     * information, but empty (not-null) fields:
     *  - Cigar String
     *  - Read Bases
     *  - Base Qualities
     *
     * Use this method if you want to create a new empty GATKSAMRecord based on
     * another GATKSAMRecord
     *
     * @param read a read to copy the header from
     * @return a read with no bases but safe for the GATK
     */
    public static GATKSAMRecord emptyRead(GATKSAMRecord read) {
        final GATKSAMRecord emptyRead = new GATKSAMRecord(read.getHeader());
        emptyRead.setReferenceIndex(read.getReferenceIndex());
        emptyRead.setAlignmentStart(0);
        emptyRead.setMappingQuality(0);
        // setting read indexing bin last
        emptyRead.setFlags(read.getFlags());
        emptyRead.setMateReferenceIndex(read.getMateReferenceIndex());
        emptyRead.setMateAlignmentStart(read.getMateAlignmentStart());
        emptyRead.setInferredInsertSize(read.getInferredInsertSize());

        emptyRead.setCigarString("");
        emptyRead.setReadBases(new byte[0]);
        emptyRead.setBaseQualities(new byte[0]);

        SAMReadGroupRecord samRG = read.getReadGroup();
        emptyRead.clearAttributes();
        if (samRG != null) {
            GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(samRG);
            emptyRead.setReadGroup(rg);
        }

        GATKBin.setReadIndexingBin(emptyRead, 0);

        return emptyRead;
    }

    /**
     * Creates a new GATKSAMRecord with the source read's header, read group and mate
     * information, but with the following fields set to user-supplied values:
     *  - Read Bases
     *  - Base Qualities
     *  - Base Insertion Qualities
     *  - Base Deletion Qualities
     *
     *  Cigar string is empty (not-null)
     *
     * Use this method if you want to create a new GATKSAMRecord based on
     * another GATKSAMRecord, but with modified bases and qualities
     *
     * @param read a read to copy the header from
     * @param readBases an array containing the new bases you wish use in place of the originals
     * @param baseQualities an array containing the new base qualities you wish use in place of the originals
     * @param baseInsertionQualities an array containing the new base insertion qaulities
     * @param baseDeletionQualities an array containing the new base deletion qualities
     * @return a read with modified bases and qualities, safe for the GATK
     */
    public static GATKSAMRecord createQualityModifiedRead(final GATKSAMRecord read,
                                                          final byte[] readBases,
                                                          final byte[] baseQualities,
                                                          final byte[] baseInsertionQualities,
                                                          final byte[] baseDeletionQualities) {
        if ( baseQualities.length != readBases.length || baseInsertionQualities.length != readBases.length || baseDeletionQualities.length != readBases.length )
            throw new IllegalArgumentException("Read bases and read quality arrays aren't the same size: Bases:" + readBases.length
                                                + " vs Base Q's:" + baseQualities.length
                                                + " vs Insert Q's:" + baseInsertionQualities.length
                                                + " vs Delete Q's:" + baseDeletionQualities.length);

        final GATKSAMRecord processedRead = GATKSAMRecord.emptyRead(read);
        processedRead.setReadBases(readBases);
        processedRead.setBaseQualities(baseQualities, EventType.BASE_SUBSTITUTION);
        processedRead.setBaseQualities(baseInsertionQualities, EventType.BASE_INSERTION);
        processedRead.setBaseQualities(baseDeletionQualities, EventType.BASE_DELETION);

        return processedRead;
    }

    /**
     * Shallow copy of everything, except for the attribute list and the temporary attributes. 
     * A new list of the attributes is created for both, but the attributes themselves are copied by reference.  
     * This should be safe because callers should never modify a mutable value returned by any of the get() methods anyway.
     * 
     * @return a shallow copy of the GATKSAMRecord
     */
    @Override
    public Object clone() {
        try {
            final GATKSAMRecord clone = (GATKSAMRecord) super.clone();
            if (temporaryAttributes != null) {
                clone.temporaryAttributes = new HashMap<>();
                for (Object attribute : temporaryAttributes.keySet())
                    clone.setTemporaryAttribute(attribute, temporaryAttributes.get(attribute));
            }
            return clone;
        } catch (final CloneNotSupportedException e) {
            throw new RuntimeException( e );
        }
    }

    /**
     * A caching version of ReadUtils.getAdaptorBoundary()
     *
     * see #ReadUtils.getAdaptorBoundary(SAMRecord) for more information about the meaning of this function
     *
     * WARNING -- this function caches a value depending on the inferred insert size and alignment starts
     * and stops of this read and its mate.  Changing these values in any way will invalidate the cached value.
     * However, we do not monitor those setter functions, so modifying a GATKSAMRecord in any way may
     * result in stale cached values.
     *
     * @return the result of calling ReadUtils.getAdaptorBoundary on this read
     */
    @Ensures("result == ReadUtils.getAdaptorBoundary(this)")
    public int getAdaptorBoundary() {
        if ( adapterBoundary == null )
            adapterBoundary = ReadUtils.getAdaptorBoundary(this);
        return adapterBoundary;
    }
}
