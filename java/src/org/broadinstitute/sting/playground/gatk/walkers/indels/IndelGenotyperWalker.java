package org.broadinstitute.sting.playground.gatk.walkers.indels;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.playground.utils.CircularArray;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;


public class IndelGenotyperWalker extends ReadWalker<Integer,Integer> {
	@Argument(fullName="bed", shortName="bed", doc="BED output file name", required=true)
	public java.io.File bed_file;
	@Argument(fullName="somatic", shortName="somatic", 
			doc="Perform somatic calls; two input alignment files must be specified", required=false)
	public boolean call_somatic = false;
	 
	private static int WINDOW_SIZE = 200;
	private RunningCoverage coverage;
	private RunningCoverage normal_coverage; // when performing somatic calls, we will be using this one for normal, and 'coverage' for tumor
	private int currentContigIndex = -1;
	private String refName = null;
	private java.io.Writer output = null;

	private Set<String> normal_samples = new HashSet<String>();
	private Set<String> tumor_samples = new HashSet<String>();
	
	@Override
	public void initialize() {
		coverage = new RunningCoverage(0,WINDOW_SIZE);
		
		int nSams = getToolkit().getArguments().samFiles.size(); 
		
		
		if ( call_somatic ) {
			if ( nSams != 2 ) {
				System.out.println("In --somatic mode two input bam files must be specified (normal/tumor)");
				System.exit(1);
			}
			normal_coverage = new RunningCoverage(0,WINDOW_SIZE);
			
			// this is an ugly hack: we want to be able to tell what file (tumor/normal sample) each read came from,
			// but reads do not carry this information!
			
			SAMFileReader rn = new SAMFileReader(getToolkit().getArguments().samFiles.get(0));
			for ( SAMReadGroupRecord rec : rn.getFileHeader().getReadGroups() ) {
				normal_samples.add(rec.getSample());
			}
			rn.close();
			rn = new SAMFileReader(getToolkit().getArguments().samFiles.get(1));
			for ( SAMReadGroupRecord rec : rn.getFileHeader().getReadGroups() ) {
				tumor_samples.add(rec.getSample());
			}
			rn.close();
			
		} else {
			if ( nSams != 1 ) System.out.println("WARNING: multiple input files specified. \n"+
					"WARNING: Without --somatic option they will be merged and processed as a single sample");
		}
		try {
			output = new java.io.FileWriter(bed_file);
		} catch (IOException e) {
			throw new StingException("Failed to open file for writing BED output");
		}
	}
	
	@Override
	public Integer map(char[] ref, SAMRecord read) {
		
		if ( read.getReadUnmappedFlag() && read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX &&
				read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ) {
			System.out.println("I think I reached unmapped reads at the end of the file(s) and I am done...");
			return -1;
		}
		
		if ( read.getReadUnmappedFlag() ||
			 read.getDuplicateReadFlag() ||
			 read.getNotPrimaryAlignmentFlag() ||
			 read.getMappingQuality() == 0 ) return 0; // we do not need those reads!
		
		
		if ( read.getReferenceIndex() != currentContigIndex ) { 
			// we just jumped onto a new contig
			
			if ( read.getReferenceIndex() < currentContigIndex ) // paranoidal 
				throw new StingException("Read "+read.getReadName()+": contig is out of order");
			
			if ( call_somatic) emit_somatic(1000000000); // print remaining indels from the previous contig (if any);
			else emit(1000000000);
			currentContigIndex = read.getReferenceIndex();
			refName = new String(read.getReferenceName());
			coverage.clear(); // reset coverage window; this will also set reference position to 0
			if ( call_somatic) normal_coverage.clear();
		}
		
		if ( read.getAlignmentStart() < coverage.getStart() ) {
			// should never happen
			throw new StingException("Read "+read.getReadName()+": out of order on the contig");
		}
		
		// reads are sorted; we are not going to see any more coverage or new indels prior
		// to current read's start position!
		
		if ( call_somatic ) emit_somatic( read.getAlignmentStart() );
		else emit( read.getAlignmentStart() );
		
		
		if ( read.getAlignmentEnd() > coverage.getStop()) {
			// should never happen
			throw new StingException("Read "+read.getReadName()+": out of coverage window bounds.\n"+
					"Read length="+read.getReadLength()+"; cigar="+read.getCigarString()+"; start="+
					read.getAlignmentStart()+"; end="+read.getAlignmentEnd()+"; window start="+coverage.getStart()+
					"; window end="+coverage.getStop());
		}
		
		if ( call_somatic ) {
			String rg = (String)read.getAttribute("RG");
			
			String sample = read.getHeader().getReadGroup(rg).getSample();
			
			if ( normal_samples.contains(sample) ) {
				normal_coverage.add(read,ref);
			} else if ( tumor_samples.contains(sample) ) {
				coverage.add(read,ref);
			} else {
				throw new StingException("Unrecognized sample: "+sample);
			}
		} else { 
			coverage.add(read, ref);
		}
		
		return 1;
	}
	
	/** Output indel calls up to the specified position and shift the coverage array(s): after this method is executed
	 * first elements of the coverage arrays map onto 'position'
	 * 
	 * @param position
	 */
	private void emit(long position) {
		
		for ( long pos = coverage.getStart() ; pos < Math.min(position,coverage.getStop()+1) ; pos++ ) {
			
			List<IndelVariant> variants = coverage.indelsAt(pos);
			if ( variants.size() == 0 ) continue; // no indels

			int cov = coverage.coverageAt(pos);
			
			if ( cov < 6 ) continue; // low coverage
			
			int total_variant_count = 0;
			int max_variant_count = 0;
			String indelString = null;
			int event_length = 0; // length of the event on the reference
			
			for ( IndelVariant var : variants ) {
				int cnt = var.getCount();
				total_variant_count +=cnt;
				if ( cnt > max_variant_count ) {
					max_variant_count = cnt;
					indelString = var.getBases();
					event_length = var.lengthOnRef();
				}
			}	
			if ( (double)total_variant_count > 0.3*cov && (double) max_variant_count > 0.7*total_variant_count ) { 		
			
				try {
					output.write(refName+"\t"+(pos-1)+"\t"+(event_length > 0 ? pos-1+event_length : pos-1)+
							"\t"+(event_length>0? "-":"+")+indelString +":"+total_variant_count+"/"+cov+"\n");
				} catch (IOException e) {
					System.out.println(e.getMessage());
					e.printStackTrace();
					throw new StingException("Error encountered while writing into output BED file");
				}
			}
//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
		}
		
		coverage.shift((int)(position - coverage.getStart() ) );
	}

	/** Output somatic indel calls up to the specified position and shift the coverage array(s): after this method is executed
	 * first elements of the coverage arrays map onto 'position'
	 * 
	 * @param position
	 */
	private void emit_somatic(long position) {
		
		for ( long pos = coverage.getStart() ; pos < Math.min(position,coverage.getStop()+1) ; pos++ ) {
			
			List<IndelVariant> tumor_variants = coverage.indelsAt(pos);
			List<IndelVariant> normal_variants = normal_coverage.indelsAt(pos);

			if ( tumor_variants.size() == 0 ) continue; // no indels in tumor

					
			int tumor_cov = coverage.coverageAt(pos);
			int normal_cov = normal_coverage.coverageAt(pos);
			
			if ( tumor_cov < 6 ) continue; // low coverage
			if ( normal_cov < 4 ) continue; // low coverage
			
			int total_variant_count_tumor = 0;
			int max_variant_count_tumor = 0;
			String indelStringTumor = null;
			int event_length_tumor = 0; // length of the event on the reference
			
			for ( IndelVariant var : tumor_variants ) {
				int cnt = var.getCount();
				total_variant_count_tumor +=cnt;
				if ( cnt > max_variant_count_tumor ) {
					max_variant_count_tumor = cnt;
					indelStringTumor = var.getBases();
					event_length_tumor = var.lengthOnRef();
				}
			}	
			
			if ( (double)total_variant_count_tumor > 0.3*tumor_cov && (double) max_variant_count_tumor > 0.7*total_variant_count_tumor ) {
				
				System.out.println(refName+"\t"+(pos-1)+"\t"+(event_length_tumor > 0 ? pos-1+event_length_tumor : pos-1)+
							"\t"+(event_length_tumor >0? "-":"+")+indelStringTumor +":"+total_variant_count_tumor+"/"+tumor_cov);
				if ( normal_variants.size() == 0 ) { 		
			
					try {
						output.write(refName+"\t"+(pos-1)+"\t"+(event_length_tumor > 0 ? pos-1+event_length_tumor : pos-1)+
								"\t"+(event_length_tumor >0? "-":"+")+indelStringTumor +":"+total_variant_count_tumor+"/"+tumor_cov+"\n");
					} catch (IOException e) {
						System.out.println(e.getMessage());
						e.printStackTrace();
						throw new StingException("Error encountered while writing into output BED file");
					}
					System.out.println("INDICATION FOR SOMATIC\n---------------------------\n");
				} else {
					System.out.println("INDICATION FOR GERMLINE\n---------------------------\n");
				}
			}
//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
		}
		
		coverage.shift((int)(position - coverage.getStart() ) );
		normal_coverage.shift((int)(position - normal_coverage.getStart() ) );
	}

	
	@Override
	public void onTraversalDone(Integer result) {
		emit(1000000000); // emit everything we might have left
		try {
			output.close();
		} catch (IOException e) {
			System.out.println("Failed to close output BED file gracefully, data may be lost");
			e.printStackTrace();
		}
		super.onTraversalDone(result);
	}
	
	@Override
	public Integer reduce(Integer value, Integer sum) {
		if ( value == -1 ) {
			onTraversalDone(sum);
			System.exit(1);
		}
		sum += value;
		return sum;
	}

	@Override
	public Integer reduceInit() {
		
		return new Integer(0);
	}
	
	
	static class IndelVariant {
		public static enum Type { I, D};
		private String bases;
		private Type type;
		private int count;
		
		public IndelVariant(Type type, String bases) {
			this.type = type;
			this.bases = bases;
			this.count = 1;
		}
		
		public void increment(int i) {
			count += i;
		}
		
		/** Returns length of the event on the reference (number of deleted bases
		 * for deletions, -1 for insertions. 
		 * @return
		 */
		public int lengthOnRef() { 
			if ( type == Type.D ) return bases.length();
			else return -1;
		}
		
		public void increment() { count+=1; }
		
		public int getCount() { return count; }
		
		public String getBases() { return bases; }
		
		public Type getType() { return type; }
		
		@Override
		public boolean equals(Object o) {
			if ( ! ( o instanceof IndelVariant ) ) return false;
			IndelVariant that = (IndelVariant)o;
			return ( this.type == that.type && this.bases.equals(that.bases) );
		}
		
		public boolean equals(Type type, String bases) {
			return ( this.type == type && this.bases.equals(bases) );
		}
	}
	
	static class RunningCoverage {
		private long start; // we keep coverage starting at this position on the reference 
		
		private CircularArray.Int coverageWindow;
		private CircularArray< List< IndelVariant > > indels;
		
		private static List<IndelVariant> emptyIndelList;
		
		static {
			emptyIndelList = new ArrayList<IndelVariant>();
		}
		
		public RunningCoverage(long start, int length) {
			this.start = start;
			coverageWindow = new CircularArray.Int(length);
			indels = new CircularArray< List<IndelVariant> >(length);
		}
		
		/** Returns 1-based reference start position of the interval this object keeps coverage for.
		 *  
		 * @return
		 */
		public long getStart() { return start; }
		
		/** Returns 1-based reference stop position (inclusive) of the interval this object keeps coverage for.
		 * 
		 * @return
		 */
		public long getStop() { return start + coverageWindow.length() - 1; }
		
		/** Returns the number of reads spanning over the specified reference position
		 * (regardless of whether they have a base or indel at that specific location)
		 * @param refPos position on the reference; must be within the bounds of the window, 
		 * otherwise IndexOutOfBoundsException will be thrown
		 */
		public int coverageAt(final long refPos) {
			return coverageWindow.get( (int)( refPos - start ) );
		}
		
		public List<IndelVariant> indelsAt( final long refPos ) {
			List<IndelVariant> l = indels.get((int)( refPos - start ));
			if ( l == null ) return emptyIndelList;
			else return l;
		}
		
		/** Increments coverage in the currently held window for every position covered by the 
		 * specified read; we count the hole span of read getAlignmentStart()-getAlignmentEnd() here,
		 * regardless of whether there are indels in the middle.Read must be completely within the current
		 * window, or an exception will be thrown.
		 * @param r
		 */
		public void add(SAMRecord r, char [] ref) {
			final long rStart = r.getAlignmentStart();
			final long rStop = r.getAlignmentEnd();

	        int localStart = (int)( rStart - start ); // start of the alignment wrt start of the current window
			
			try {
				for ( int k = localStart; k <= (int)(rStop-start) ; k++ ) coverageWindow.increment(k, 1);
			} catch ( IndexOutOfBoundsException e) { // replace the message and re-throw:
				throw new IndexOutOfBoundsException("Current coverage window: "+getStart()+"-"+getStop()+
						"; illegal attempt to add read spanning "+rStart+"-"+rStop);				
			}
			
			// now let's extract indels:
			
		   	Cigar c = r.getCigar();
	    	final int nCigarElems = c.numCigarElements();

	    	// if read has no indels, there is nothing to do
	        if ( c.numCigarElements() <= 1 ) return ; 
	    	
	        int posOnRead = 0;
	        int posOnRef = 0; // the chunk of reference ref[] that we have access to is aligned with the read:
	        	              // its start on the actual full reference contig is r.getAlignmentStart()
	        
	        for ( int i = 0 ; i < nCigarElems ; i++ ) {

	            final CigarElement ce = c.getCigarElement(i);
	            IndelVariant.Type type = null;
	            String bases = null;
	            
	            int indelPosition = 0; // indel position in our coverage window (i.e. relative to getStart()).
	            	                   // note that here we assign indels to the first deleted base or to the first
	            	                   // base after insertion
	            
	            switch(ce.getOperator()) {
	            case I:
	                    type = IndelVariant.Type.I; 
	                    bases = r.getReadString().substring(posOnRead,posOnRead+ce.getLength()); 
	                    indelPosition = localStart + posOnRef ;
	                    // will increment position on the read below, there's no 'break' statement yet...
	            case H:
	            case S:
	            		// here we also skip hard and soft-clipped bases on th eread; according to SAM format specification, 
               			// alignment start position on the reference points to where the actually aligned 
    					// (not clipped) bases go, so we do not need to increment reference position here	                    
	            		posOnRead += ce.getLength();
	                    break;
	            case D: 
	            		type = IndelVariant.Type.D;
	            		bases = new String( ref, posOnRef, ce.getLength() );
	                    indelPosition = localStart + posOnRef ;
	                    posOnRef += ce.getLength();
	                    break;
	            case M: 
	            		posOnRef += ce.getLength(); 
	            		posOnRead += ce.getLength(); 
	            		break; // advance along the gapless block in the alignment
	            default :
	                throw new IllegalArgumentException("Unexpected operator in cigar string: "+ce.getOperator());
	            }

	            if ( type == null ) continue; // element was not an indel, go grab next element...
	            
	            // we got an indel if we are here...
	            if ( i == 0 ) logger.warn("Indel at the start of the read "+r.getReadName());
	            if ( i == nCigarElems - 1) logger.warn("Indel at the end of the read "+r.getReadName());

	            updateCount(indelPosition, type, bases);
	        }
	            
			
		}
		
		/** Convenience shortcut method. Checks if indel of specified type and with specified bases is already recorded
		 * for position <code>pos</code> (relative to start of the window getStart()). If such indel is found, the counter
		 * is increased; if it is not found, a new indel (with count = 1, obviously) will be added at that position. If indel array
		 * still had null at the specified position, this method will instantiate new list of indels for this position
		 * transparently.
		 * 
		 * @param pos
		 * @param type
		 * @param bases
		 */
		private void updateCount(int pos, IndelVariant.Type type, String bases) {
            List<IndelVariant> indelsAtSite = indels.get(pos);
            if ( indelsAtSite == null ) {
            	indelsAtSite = new ArrayList<IndelVariant>();
            	indels.set(pos, indelsAtSite);
            }
            
            boolean found = false;
            for ( IndelVariant v : indelsAtSite ) {
            	if ( ! v.equals(type, bases) ) continue;

            	v.increment();
            	found = true;
            	break;
            }
            
            if ( ! found ) indelsAtSite.add(new IndelVariant(type, bases));
			
		}
		
		/** Resets reference start position to 0 and sets all coverage counts in the window to 0.
		 * 
		 */
		public void clear() { 
			start = 0; 
			coverageWindow.clear();
			indels.clear();
		}
		
		/** Shifts current window to the right along the reference contig by the specified number of bases.
		 * Coverage counts computed earlier for the positions that remain in scope will be preserved.
		 * @param offset
		 */
		public void shift(int offset) {
			start += offset;
			coverageWindow.shiftData(offset);
			indels.shiftData(offset);
		}
	}

}
