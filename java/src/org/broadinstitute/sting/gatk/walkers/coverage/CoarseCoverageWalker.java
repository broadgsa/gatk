package org.broadinstitute.sting.gatk.walkers.coverage;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;

/**
 * Computes the coverage per every <granularity> bases on the reference, or on the subset of the reference
 * specified by the intervals provided. Moving to the next contig on the reference will always restart the
 * count anew, even if the count of bases in the last chunk on the previous contig did not reach specified <granularity>.
 */
public class CoarseCoverageWalker extends ReadWalker<Integer,Integer> {
	
    @Argument(fullName="granularity", shortName="G", doc="Granularity", required=true)
    public Integer N;
	
    private int chunkStart = 1; // start of the current chunk we are counting reads for
    private int contig = 0; // current contig we are on
    private int count = 0; // number of reads overlapping with the current chunk
    private static String zeroString = "0";
    
    @Override 
    public void initialize() {
    	chunkStart = 1;
    	contig = 0;
    	count = 0;
    }
    
	@Override
	public Integer map(char[] ref, SAMRecord read) {
		
		if ( read.getReadUnmappedFlag() ||
			 read.getDuplicateReadFlag() ||
			 read.getNotPrimaryAlignmentFlag() ||
			 read.getMappingQuality() == 0 ) return 0;
		
		if ( read.getReferenceIndex() != contig ) {
			// we jumped onto another contig
			out.printf("%d%n", count); // print old count
			count = 0;
			
			// if we skipped one or more contigs completely, make sure we print 0 counts over all of them:
			for ( contig++ ; contig < read.getReferenceIndex() ; contig++) {
				int contigSize = read.getHeader().getSequence(contig).getSequenceLength();
				for ( int k = 1 ; k < contigSize ;  k+=N ) out.println(zeroString);
			}
			// by now we scrolled to the right contig
			
			chunkStart = 1; // reset chunk start
		}
		
		// if our read is past the boundary of the current chunk, print old count(s)
		// (for the current chunk and all chunks we may have skipped altogether) and reinitialize:
		while ( chunkStart+N < read.getAlignmentStart() ) {
			out.printf("%d%n", count); // print old count
			count = 0;
			chunkStart += N;
		}
		count++;
		return 1;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		return value+sum;
	}

	@Override
	public Integer reduceInit() {
		return 0;
	}
	
	@Override 
	public void onTraversalDone(Integer result) {
		out.printf("%d%n", count); // print count from the last chunk
		super.onTraversalDone(result);
		
	}
	

}
