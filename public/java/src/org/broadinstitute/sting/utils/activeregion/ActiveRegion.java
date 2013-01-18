/*
* Copyright (c) 2012 The Broad Institute
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
        extendedLoc = genomeLocParser.createGenomeLocOnContig(activeRegionLoc.getContig(), activeRegionLoc.getStart() - extension, activeRegionLoc.getStop() + extension);
        fullExtentReferenceLoc = extendedLoc;
    }

    @Override
    public String toString() {
        return "ActiveRegion "  + activeRegionLoc.toString() + " active?=" + isActive + " nReads=" + reads.size() + " ";
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
}