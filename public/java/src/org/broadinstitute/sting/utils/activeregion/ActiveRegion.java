package org.broadinstitute.sting.utils.activeregion;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 1/4/12
 */

public class ActiveRegion implements HasGenomeLocation {

    private final ArrayList<ActiveRead> reads = new ArrayList<ActiveRead>();
    private byte[] reference = null;
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

    // add each read to the bin and extend the reference genome activeRegionLoc if needed
    public void add( final GATKSAMRecord read, final boolean isPrimaryRegion  ) {
        fullExtentReferenceLoc = fullExtentReferenceLoc.union( genomeLocParser.createGenomeLoc( read ) );
        reads.add( new ActiveRead(read, isPrimaryRegion) );
    }

    public ArrayList<ActiveRead> getReads() { return reads; }

    public byte[] getReference( final IndexedFastaSequenceFile referenceReader ) {
        // set up the reference if we haven't done so yet
        if ( reference == null ) {
            reference = referenceReader.getSubsequenceAt(fullExtentReferenceLoc.getContig(), fullExtentReferenceLoc.getStart(), fullExtentReferenceLoc.getStop()).getBases();
        }

        return reference;
    }

    @Override
    public GenomeLoc getLocation() { return activeRegionLoc; }

    public GenomeLoc getExtendedLoc() { return extendedLoc; }
    public GenomeLoc getReferenceLoc() { return fullExtentReferenceLoc; }

    public int getExtension() { return extension; }
    public int size() { return reads.size(); }
}