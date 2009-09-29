package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReadBackedPileup extends BasicPileup {
    GenomeLoc loc;
    char ref;
    List<SAMRecord> reads;
    List<Integer> offsets;

    public ReadBackedPileup(char ref, AlignmentContext context ) {
        this(context.getLocation(), ref, context.getReads(), context.getOffsets());
    }

    public ReadBackedPileup(GenomeLoc loc, char ref, List<SAMRecord> reads, List<Integer> offsets ) {
        assert reads != null;
        assert offsets != null;
        assert reads.size() == offsets.size();
        
        this.loc = loc;
        this.ref = ref;
        this.reads = reads;
        this.offsets = offsets;
    }

    public int size()                 { return reads.size(); }
    public List<SAMRecord> getReads() { return reads; }
    public List<Integer> getOffsets() { return offsets; }

    public GenomeLoc getLocation() {
        return loc;
    }

    public char getRef() {
        return ref;
    }

    public String getBases() {
        return basePileupAsString(reads, offsets, includeDeletions);
    }

    public String getBasesWithStrand() {
        return baseWithStrandPileupAsString(reads, offsets, includeDeletions);
    }

    public String getQuals() {
        return qualPileupAsString(reads, offsets);
    }

    public String getQualsAsInts() {
        //System.out.printf("getQualsAsInts");        
        return Utils.join(",", qualPileup(reads, offsets));
    }

    public String getMappingQualsAsInts() {
        return Utils.join(",", mappingQualPileup(reads));
    }

    public String getMappingQuals() {
        return mappingQualPileupAsString(reads);
    }

    public String getSecondaryBasePileup() {
        return secondaryBasePileupAsString(reads, offsets);
    }

    public String getSecondaryQualPileup() {
        return secondaryQualPileupAsString(reads, offsets);
    }

    public String getBasePileupAsCountsString() {
        String bases = basePileupAsString(reads, offsets, includeDeletions);

	int[] counts = new int[4];
	for (int i = 0; i < reads.size(); i++)
	    {
		// skip deletion sites
		if ( offsets.get(i) == -1 )
		    continue;
		char base = Character.toUpperCase((char)(reads.get(i).getReadBases()[offsets.get(i)]));
		if (BaseUtils.simpleBaseToBaseIndex(base) == -1) { continue; }
		counts[BaseUtils.simpleBaseToBaseIndex(base)]++; 
	    }
	return String.format("A[%d] C[%d] G[%d] T[%d]",
			     counts[0],
			     counts[1],
			     counts[2],
			     counts[3]);
    }

    public String getProbDistPileup() {
        return probDistPileupAsString(reads, offsets);
    }

    public String toString() {
        return getPileupString(true);
    }

    public String getPileupString(boolean qualsAsInts)
    {
        // In the pileup format, each line represents a genomic position, consisting of chromosome name,
        // coordinate, reference base, read bases, read qualities and alignment mapping qualities.

        //return String.format("%s %s %s %s", getLocation(), getRef(), getBases(), getQuals());
        return String.format("%s %s %s %s %s %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                getRef(),                                               // reference base
                getBases(),
                qualsAsInts ? getQualsAsInts() : getQuals(),
                qualsAsInts ? getMappingQualsAsInts() : getMappingQuals() );
    }

    public String getPileupWithStrandString(boolean qualsAsInts) {
        return String.format("%s %s %s %s %s %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                getRef(),                                               // reference base
                getBasesWithStrand(),
                qualsAsInts ? getQualsAsInts() : getQuals(),
                qualsAsInts ? getMappingQualsAsInts() : getMappingQuals() );
    }
}
