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

package org.broadinstitute.gatk.utils.activeregion;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.HasGenomeLocation;
import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.util.*;

/**
 * Represents a single active region created by the Active Region Traversal for processing
 *
 * An active region is a single contiguous span of bases on the genome that should be operated
 * on as a single unit for the active region traversal.  The action may contains a list of
 * reads that overlap the region (may because there may be no reads in the region).  The region
 * is tagged as being either active or inactive, depending on the probabilities provided by
 * the isActiveProb results from the ART walker.  Each region carries with it the
 * exact span of the region (bases which are the core of the isActiveProbs from the walker) as
 * well as an extended size, that includes the ART walker's extension size.  Reads in the region
 * provided by ART include all reads overlapping the extended span, not the raw span.
 *
 * User: rpoplin
 * Date: 1/4/12
 */
@Invariant({
        "extension >= 0",
        "activeRegionLoc != null",
        "genomeLocParser != null",
        "spanIncludingReads != null",
        "extendedLoc != null"
})
public class ActiveRegion implements HasGenomeLocation {
    /**
     * The reads included in this active region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    private final List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();

    /**
     * An ordered list (by genomic coordinate) of the ActivityProfileStates that went
     * into this active region.  May be empty, which says that no supporting states were
     * provided when this region was created.
     */
    private final List<ActivityProfileState> supportingStates;

    /**
     * The raw span of this active region, not including the active region extension
     */
    private final GenomeLoc activeRegionLoc;

    /**
     * The span of this active region on the genome, including the active region extension
     */
    private final GenomeLoc extendedLoc;

    /**
     * The extension, in bp, of this active region.
     */
    private final int extension;

    /**
     * A genomeLocParser so we can create genomeLocs
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    private final boolean isActive;

    /**
     * The span of this active region, including the bp covered by all reads in this
     * region.  This union of extensionLoc and the loc of all reads in this region.
     *
     * Must be at least as large as extendedLoc, but may be larger when reads
     * partially overlap this region.
     */
    private GenomeLoc spanIncludingReads;


    /**
     * Indicates whether the active region has been finalized
     */
    private boolean hasBeenFinalized;

    /**
     * Create a new ActiveRegion containing no reads
     *
     * @param activeRegionLoc the span of this active region
     * @param supportingStates the states that went into creating this region, or null / empty if none are available.
     *                         If not empty, must have exactly one state for each bp in activeRegionLoc
     * @param isActive indicates whether this is an active region, or an inactve one
     * @param genomeLocParser a non-null parser to let us create new genome locs
     * @param extension the active region extension to use for this active region
     */
    public ActiveRegion( final GenomeLoc activeRegionLoc, final List<ActivityProfileState> supportingStates, final boolean isActive, final GenomeLocParser genomeLocParser, final int extension ) {
        if ( activeRegionLoc == null ) throw new IllegalArgumentException("activeRegionLoc cannot be null");
        if ( activeRegionLoc.size() == 0 ) throw new IllegalArgumentException("Active region cannot be of zero size, but got " + activeRegionLoc);
        if ( genomeLocParser == null ) throw new IllegalArgumentException("genomeLocParser cannot be null");
        if ( extension < 0 ) throw new IllegalArgumentException("extension cannot be < 0 but got " + extension);

        this.activeRegionLoc = activeRegionLoc;
        this.supportingStates = supportingStates == null ? Collections.<ActivityProfileState>emptyList() : new ArrayList<ActivityProfileState>(supportingStates);
        this.isActive = isActive;
        this.genomeLocParser = genomeLocParser;
        this.extension = extension;
        this.extendedLoc = genomeLocParser.createGenomeLocOnContig(activeRegionLoc.getContig(), activeRegionLoc.getStart() - extension, activeRegionLoc.getStop() + extension);
        this.spanIncludingReads = extendedLoc;

        if ( ! this.supportingStates.isEmpty() ) {
            if ( this.supportingStates.size() != activeRegionLoc.size() )
                throw new IllegalArgumentException("Supporting states wasn't empty but it doesn't have exactly one state per bp in the active region: states " + this.supportingStates.size() + " vs. bp in region = " + activeRegionLoc.size());
            GenomeLoc lastStateLoc = null;
            for ( final ActivityProfileState state : this.supportingStates ) {
                if ( lastStateLoc != null ) {
                    if ( state.getLoc().getStart() != lastStateLoc.getStart() + 1 || state.getLoc().getContigIndex() != lastStateLoc.getContigIndex())
                        throw new IllegalArgumentException("Supporting state has an invalid sequence: last state was " + lastStateLoc + " but next state was " + state);
                }
                lastStateLoc = state.getLoc();
            }
        }
    }

    /**
     * Simple interface to create an active region that isActive without any profile state
     */
    public ActiveRegion( final GenomeLoc activeRegionLoc, final GenomeLocParser genomeLocParser, final int extension ) {
        this(activeRegionLoc, Collections.<ActivityProfileState>emptyList(), true, genomeLocParser, extension);
    }

    @Override
    public String toString() {
        return "ActiveRegion "  + activeRegionLoc.toString() + " active?=" + isActive() + " nReads=" + reads.size();
    }

    /**
     * See #getActiveRegionReference but with padding == 0
     */
    public byte[] getActiveRegionReference( final IndexedFastaSequenceFile referenceReader ) {
        return getActiveRegionReference(referenceReader, 0);
    }

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this active region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    @Ensures("result != null")
    public byte[] getActiveRegionReference( final IndexedFastaSequenceFile referenceReader, final int padding ) {
        return getReference(referenceReader, padding, extendedLoc);
    }

    /**
     * See #getActiveRegionReference but using the span including regions not the extended span
     */
    public byte[] getFullReference( final IndexedFastaSequenceFile referenceReader ) {
        return getFullReference(referenceReader, 0);
    }

    /**
     * See #getActiveRegionReference but using the span including regions not the extended span
     */
    public byte[] getFullReference( final IndexedFastaSequenceFile referenceReader, final int padding ) {
        return getReference(referenceReader, padding, spanIncludingReads);
    }

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this active region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @param genomeLoc a non-null genome loc indicating the base span of the bp we'd like to get the reference for
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    @Ensures("result != null")
    public byte[] getReference( final IndexedFastaSequenceFile referenceReader, final int padding, final GenomeLoc genomeLoc ) {
        if ( referenceReader == null ) throw new IllegalArgumentException("referenceReader cannot be null");
        if ( padding < 0 ) throw new IllegalArgumentException("padding must be a positive integer but got " + padding);
        if ( genomeLoc == null ) throw new IllegalArgumentException("genomeLoc cannot be null");
        if ( genomeLoc.size() == 0 ) throw new IllegalArgumentException("GenomeLoc must have size > 0 but got " + genomeLoc);

        final byte[] reference =  referenceReader.getSubsequenceAt( genomeLoc.getContig(),
                Math.max(1, genomeLoc.getStart() - padding),
                Math.min(referenceReader.getSequenceDictionary().getSequence(genomeLoc.getContig()).getSequenceLength(), genomeLoc.getStop() + padding) ).getBases();

        return reference;
    }

    /**
     * Get the raw span of this active region (excluding the extension)
     * @return a non-null genome loc
     */
    @Override
    @Ensures("result != null")
    public GenomeLoc getLocation() { return activeRegionLoc; }

    /**
     * Get the span of this active region including the extension value
     * @return a non-null GenomeLoc
     */
    @Ensures("result != null")
    public GenomeLoc getExtendedLoc() { return extendedLoc; }

    /**
     * Get the span of this active region including the extension and the projects on the
     * genome of all reads in this active region.  That is, returns the bp covered by this
     * region and all reads in the region.
     * @return a non-null genome loc
     */
    @Ensures("result != null")
    public GenomeLoc getReadSpanLoc() { return spanIncludingReads; }

    /**
     * Get the active profile states that went into creating this region, if possible
     * @return an unmodifiable list of states that led to the creation of this region, or an empty
     *         list if none were provided
     */
    @Ensures("result != null")
    public List<ActivityProfileState> getSupportingStates() {
        return Collections.unmodifiableList(supportingStates);
    }

    /**
     * Get the active region extension applied to this region
     *
     * The extension is >= 0 bp in size, and indicates how much padding this art walker wanted for its regions
     *
     * @return the size in bp of the region extension
     */
    @Ensures("result >= 0")
    public int getExtension() { return extension; }

    /**
     * Get an unmodifiable list of reads currently in this active region.
     *
     * The reads are sorted by their coordinate position
     *
     * @return an unmodifiable list of reads in this active region
     */
    @Ensures("result != null")
    public List<GATKSAMRecord> getReads() {
        return Collections.unmodifiableList(reads);
    }

    /**
     * Get the number of reads currently in this active region
     * @return an integer >= 0
     */
    @Ensures("result >= 0")
    public int size() { return reads.size(); }

    /**
     * Add read to this active region
     *
     * Read must have alignment start >= than the last read currently in this active region.
     *
     * @throws IllegalArgumentException if read doesn't overlap the extended region of this active region
     *
     * @param read a non-null GATKSAMRecord
     */
    @Ensures("reads.size() == old(reads.size()) + 1")
    public void add( final GATKSAMRecord read ) {
        if ( read == null ) throw new IllegalArgumentException("Read cannot be null");

        final GenomeLoc readLoc = genomeLocParser.createGenomeLoc( read );
        if ( ! readOverlapsRegion(read) )
            throw new IllegalArgumentException("Read location " + readLoc + " doesn't overlap with active region extended span " + extendedLoc);

        spanIncludingReads = spanIncludingReads.union( readLoc );

        if ( ! reads.isEmpty() ) {
            final GATKSAMRecord lastRead = reads.get(size() - 1);
            if ( ! lastRead.getReferenceIndex().equals(read.getReferenceIndex()) )
                throw new IllegalArgumentException("Attempting to add a read to ActiveRegion not on the same contig as other reads: lastRead " + lastRead + " attempting to add " + read);

            if ( read.getAlignmentStart() < lastRead.getAlignmentStart() )
                throw new IllegalArgumentException("Attempting to add a read to ActiveRegion out of order w.r.t. other reads: lastRead " + lastRead + " at " + lastRead.getAlignmentStart() + " attempting to add " + read + " at " + read.getAlignmentStart());
        }

        reads.add( read );
    }

    /**
     * Returns true if read would overlap the extended extent of this region
     * @param read the read we want to test
     * @return true if read can be added to this region, false otherwise
     */
    public boolean readOverlapsRegion(final GATKSAMRecord read) {
        final GenomeLoc readLoc = genomeLocParser.createGenomeLoc( read );
        return readLoc.overlapsP(extendedLoc);
    }

    /**
     * Add all reads to this active region
     * @param reads a collection of reads to add to this active region
     */
    public void addAll(final Collection<GATKSAMRecord> reads) {
        if ( reads == null ) throw new IllegalArgumentException("reads cannot be null");
        for ( final GATKSAMRecord read : reads )
            add(read);
    }

    /**
     * Clear all of the reads currently in this active region
     */
    @Ensures("size() == 0")
    public void clearReads() {
        spanIncludingReads = extendedLoc;
        reads.clear();
    }

    /**
     * Remove all of the reads in readsToRemove from this active region
     * @param readsToRemove the set of reads we want to remove
     */
    public void removeAll( final Set<GATKSAMRecord> readsToRemove ) {
        final Iterator<GATKSAMRecord> it = reads.iterator();
        spanIncludingReads = extendedLoc;
        while ( it.hasNext() ) {
            final GATKSAMRecord read = it.next();
            if ( readsToRemove.contains(read) )
                it.remove();
            else
                spanIncludingReads = spanIncludingReads.union( genomeLocParser.createGenomeLoc(read) );
        }
    }

    /**
     * Is this region equal to other, excluding any reads in either region in the comparison
     * @param other the other active region we want to test
     * @return true if this region is equal, excluding any reads and derived values, to other
     */
    protected boolean equalExceptReads(final ActiveRegion other) {
        if ( activeRegionLoc.compareTo(other.activeRegionLoc) != 0 ) return false;
        if ( isActive() != other.isActive()) return false;
        if ( genomeLocParser != other.genomeLocParser ) return false;
        if ( extension != other.extension ) return false;
        if ( extendedLoc.compareTo(other.extendedLoc) != 0 ) return false;
        return true;
    }

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    public boolean isActive() {
        return isActive;
    }

    /**
     * Intersect this active region with the allowed intervals, returning a list of active regions
     * that only contain locations present in intervals
     *
     * Note that the returned list may be empty, if this active region doesn't overlap the set at all
     *
     * Note that the resulting regions are all empty, regardless of whether the current active region has reads
     *
     * @param intervals a non-null set of intervals that are allowed
     * @return an ordered list of active region where each interval is contained within intervals
     */
    @Ensures("result != null")
    protected List<ActiveRegion> splitAndTrimToIntervals(final GenomeLocSortedSet intervals) {
        final List<GenomeLoc> allOverlapping = intervals.getOverlapping(getLocation());
        final List<ActiveRegion> clippedRegions = new LinkedList<ActiveRegion>();

        for ( final GenomeLoc overlapping : allOverlapping ) {
            clippedRegions.add(trim(overlapping, extension));
        }

        return clippedRegions;
    }

    /**
     * Trim this active to just the span, producing a new active region without any reads that has only
     * the extent of newExtend intersected with the current extent
     * @param span the new extend of the active region we want
     * @param extension the extension size we want for the newly trimmed active region
     * @return a non-null, empty active region
     */
    public ActiveRegion trim(final GenomeLoc span, final int extension) {
        if ( span == null ) throw new IllegalArgumentException("Active region extent cannot be null");
        if ( extension < 0) throw new IllegalArgumentException("the extension size must be 0 or greater");
        final int extendStart = Math.max(1,span.getStart() - extension);
        final int maxStop = genomeLocParser.getContigs().getSequence(span.getContigIndex()).getSequenceLength();
        final int extendStop = Math.min(span.getStop() + extension, maxStop);
        final GenomeLoc extendedSpan = genomeLocParser.createGenomeLoc(span.getContig(), extendStart, extendStop);
        return trim(span, extendedSpan);

//TODO - Inconsiste support of substates trimming. Check lack of consistency!!!!
//        final GenomeLoc subLoc = getLocation().intersect(span);
//        final int subStart = subLoc.getStart() - getLocation().getStart();
//        final int subEnd = subStart + subLoc.size();
//        final List<ActivityProfileState> subStates = supportingStates.isEmpty() ? supportingStates : supportingStates.subList(subStart, subEnd);
//        return new ActiveRegion( subLoc, subStates, isActive, genomeLocParser, extension );

    }

    public ActiveRegion trim(final GenomeLoc span) {
        return trim(span,span);
    }

    /**
     * Trim this active to no more than the span, producing a new active region with properly trimmed reads that
     * attempts to provide the best possible representation of this active region covering the span.
     *
     * The challenge here is that span may (1) be larger than can be represented by this active region
     * + its original extension and (2) the extension must be symmetric on both sides.  This algorithm
     * therefore determines how best to represent span as a subset of the span of this
     * region with a padding value that captures as much of the span as possible.
     *
     * For example, suppose this active region is
     *
     * Active:    100-200 with extension of 50, so that the true span is 50-250
     * NewExtent: 150-225 saying that we'd ideally like to just have bases 150-225
     *
     * Here we represent the active region as a active region from 150-200 with 25 bp of padding.
     *
     * The overall constraint is that the active region can never exceed the original active region, and
     * the extension is chosen to maximize overlap with the desired region
     *
     * @param span the new extend of the active region we want
     * @return a non-null, empty active region
     */
    public ActiveRegion trim(final GenomeLoc span, final GenomeLoc extendedSpan) {
        if ( span == null ) throw new IllegalArgumentException("Active region extent cannot be null");
        if ( extendedSpan == null ) throw new IllegalArgumentException("Active region extended span cannot be null");
        if ( ! extendedSpan.containsP(span))
            throw new IllegalArgumentException("The requested extended must fully contain the requested span");

        final GenomeLoc subActive = getLocation().intersect(span);
        final int requiredOnRight = Math.max(extendedSpan.getStop() - subActive.getStop(), 0);
        final int requiredOnLeft = Math.max(subActive.getStart() - extendedSpan.getStart(), 0);
        final int requiredExtension = Math.min(Math.max(requiredOnLeft, requiredOnRight), getExtension());

        final ActiveRegion result = new ActiveRegion( subActive, Collections.<ActivityProfileState>emptyList(), isActive, genomeLocParser, requiredExtension );

        final List<GATKSAMRecord> myReads = getReads();
        final GenomeLoc resultExtendedLoc = result.getExtendedLoc();
        final int resultExtendedLocStart = resultExtendedLoc.getStart();
        final int resultExtendedLocStop = resultExtendedLoc.getStop();

        final List<GATKSAMRecord> trimmedReads = new ArrayList<>(myReads.size());
        for( final GATKSAMRecord read : myReads ) {
            final GATKSAMRecord clippedRead = ReadClipper.hardClipToRegion(read,
                    resultExtendedLocStart, resultExtendedLocStop);
            if( result.readOverlapsRegion(clippedRead) && clippedRead.getReadLength() > 0 )
                trimmedReads.add(clippedRead);
        }
        result.clearReads();
        result.addAll(ReadUtils.sortReadsByCoordinate(trimmedReads));
        return result;
    }

    public void setFinalized(final boolean value) {
        hasBeenFinalized = value;
    }

    public boolean isFinalized() {
        return hasBeenFinalized;
    }

}