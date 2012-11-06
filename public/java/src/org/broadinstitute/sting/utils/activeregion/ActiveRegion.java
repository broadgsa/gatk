package org.broadinstitute.sting.utils.activeregion;

import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 1/4/12
 */

public class ActiveRegion implements HasGenomeLocation {

    private final ArrayList<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
    private final GenomeLoc activeRegionLoc;
    private final GenomeLoc extendedLoc;
    private final int extension;
    private GenomeLoc fullExtentReferenceLoc = null;
    private final GenomeLocParser genomeLocParser;
    public final boolean isActive;

    public ActiveRegion( final GenomeLoc activeRegionLoc, final boolean isActive, final GenomeLocParser genomeLocParser, final int extension ) {
        this.activeRegionLoc = activeRegionLoc;
        this.isActive = isActive;
        this.genomeLocParser = genomeLocParser;
        this.extension = extension;
        extendedLoc = genomeLocParser.createGenomeLoc(activeRegionLoc.getContig(), activeRegionLoc.getStart() - extension, activeRegionLoc.getStop() + extension);
        fullExtentReferenceLoc = extendedLoc;
    }

    @Override
    public String toString() {
        return "ActiveRegion " + activeRegionLoc.toString();
    }

    // add each read to the bin and extend the reference genome activeRegionLoc if needed
    public void add( final GATKSAMRecord read ) {
        fullExtentReferenceLoc = fullExtentReferenceLoc.union( genomeLocParser.createGenomeLoc( read ) );
        reads.add( read );
    }
    
    public void hardClipToActiveRegion() {
        final ArrayList<GATKSAMRecord> clippedReads = ReadClipper.hardClipToRegion( reads, extendedLoc.getStart(), extendedLoc.getStop() );
        reads.clear();
        reads.addAll(clippedReads);
    }

    public ArrayList<GATKSAMRecord> getReads() { return reads; }

    @Requires("referenceReader.isUppercasingBases()")
    public byte[] getActiveRegionReference( final CachingIndexedFastaSequenceFile referenceReader ) {
        return getActiveRegionReference(referenceReader, 0);
    }

    @Requires("referenceReader.isUppercasingBases()")
    public byte[] getActiveRegionReference( final CachingIndexedFastaSequenceFile referenceReader, final int padding ) {
        return getReference( referenceReader, padding, extendedLoc );
    }

    @Requires("referenceReader.isUppercasingBases()")
    public byte[] getFullReference( final CachingIndexedFastaSequenceFile referenceReader ) {
        return getFullReference(referenceReader, 0);
    }

    @Requires("referenceReader.isUppercasingBases()")
    public byte[] getFullReference( final CachingIndexedFastaSequenceFile referenceReader, final int padding ) {
        return getReference( referenceReader, padding, fullExtentReferenceLoc );
    }

    @Requires("referenceReader.isUppercasingBases()")
    private byte[] getReference( final CachingIndexedFastaSequenceFile referenceReader, final int padding, final GenomeLoc genomeLoc ) {
        final byte[] reference =  referenceReader.getSubsequenceAt( genomeLoc.getContig(),
                Math.max(1, genomeLoc.getStart() - padding),
                Math.min(referenceReader.getSequenceDictionary().getSequence(genomeLoc.getContig()).getSequenceLength(), genomeLoc.getStop() + padding) ).getBases();
        return reference;
    }

    @Override
    public GenomeLoc getLocation() { return activeRegionLoc; }
    public GenomeLoc getExtendedLoc() { return extendedLoc; }
    public GenomeLoc getReferenceLoc() { return fullExtentReferenceLoc; }

    public int getExtension() { return extension; }
    public int size() { return reads.size(); }
    public void clearReads() { reads.clear(); }
    public void remove( final GATKSAMRecord read ) { reads.remove( read ); }
    public void removeAll( final ArrayList<GATKSAMRecord> readsToRemove ) { reads.removeAll( readsToRemove ); }

    public boolean equalExceptReads(final ActiveRegion other) {
        if ( activeRegionLoc.compareTo(other.activeRegionLoc) != 0 ) return false;
        if ( isActive != other.isActive ) return false;
        if ( genomeLocParser != other.genomeLocParser ) return false;
        if ( extension != other.extension ) return false;
        if ( extendedLoc.compareTo(other.extendedLoc) != 0 ) return false;
        return true;
    }

    /**
     * A comparator class which is used to sort ActiveRegions by their start location
     */
    /*
    public static class ActiveRegionStartLocationComparator implements Comparator<ActiveRegion> {

        public ActiveRegionStartLocationComparator() {}

        @Override
        public int compare(final ActiveRegion left, final ActiveRegion right) {
            return left.getLocation().compareTo(right.getLocation());
        }
    }
    */
}