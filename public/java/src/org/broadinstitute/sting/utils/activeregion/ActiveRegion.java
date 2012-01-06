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
    private final GenomeLoc loc;
    private GenomeLoc referenceLoc = null;
    private final GenomeLocParser genomeLocParser;
    public final boolean isActive;

    public ActiveRegion( final GenomeLoc loc, final boolean isActive, final GenomeLocParser genomeLocParser ) {
        this.loc = loc;
        this.isActive = isActive;
        this.genomeLocParser = genomeLocParser;
        referenceLoc = loc;
    }

    // add each read to the bin and extend the reference genome loc if needed
    public void add( final GATKSAMRecord read, final boolean isPrimaryRegion  ) {
        referenceLoc = referenceLoc.union( genomeLocParser.createGenomeLoc( read ) );
        reads.add( new ActiveRead(read, isPrimaryRegion) );
    }

    public ArrayList<ActiveRead> getReads() { return reads; }

    public byte[] getReference( final IndexedFastaSequenceFile referenceReader ) {
        // set up the reference if we haven't done so yet
        if ( reference == null ) {
            reference = referenceReader.getSubsequenceAt(referenceLoc.getContig(), referenceLoc.getStart(), referenceLoc.getStop()).getBases();
        }

        return reference;
    }

    public GenomeLoc getLocation() { return loc; }
    
    public GenomeLoc getReferenceLocation() { return referenceLoc; }

    public int size() { return reads.size(); }
}