package org.broadinstitute.sting.playground.gatk.walkers.indels;

import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.playground.utils.CircularArray;
import org.broadinstitute.sting.utils.StingException;


public class IndelGenotyperWalker extends ReadWalker<Integer,Integer> {

	private RunningCoverage coverage;
	private int currentContigIndex = -1;
	
	@Override
	public void initialize() {
		coverage = new RunningCoverage(0,100);
	}
	
	@Override
	public Integer map(char[] ref, SAMRecord read) {
		
		if ( read.getReadUnmappedFlag() ||
			 read.getDuplicateReadFlag() ||
			 read.getNotPrimaryAlignmentFlag() ||
			 read.getMappingQuality() == 0 ) return 0; // we do not need those reads!
		
		if ( read.getReferenceIndex() != currentContigIndex ) {
			
			if ( read.getReferenceIndex() < currentContigIndex ) // paranoidal 
				throw new StingException("Read "+read.getReadName()+": contig is out of order");
			
			currentContigIndex = read.getReferenceIndex();
			coverage.clear(); // reset coverage window; this will also set reference position to 0
		}
		
		if ( read.getAlignmentStart() < coverage.getStart() ) {
			// should never happen
			throw new StingException("Read "+read.getReadName()+": out of order on the contig");
		}
		
		// reads are sorted; we are not going to see any more coverage or new indels prior
		// to current read's start position!
		for ( long pos = coverage.getStart() ; pos < Math.min(read.getAlignmentStart(),coverage.getStop()+1) ; pos++ ) {
			List<IndelVariant> variants = coverage.indelsAt(pos);
			if ( variants.size() == 0 ) continue;
			System.out.print(read.getReferenceName()+"\t"+pos+"\t"+coverage.coverageAt(pos));
			for ( IndelVariant var : variants ) {
				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
			}
			System.out.println();
		}
		
		coverage.shift((int)(read.getAlignmentStart() - coverage.getStart() ) );
		
		if ( read.getAlignmentEnd() > coverage.getStop()) {
			// should never happen
			throw new StingException("Read "+read.getReadName()+": out of coverage window bounds");
		}
		
		coverage.add(read,ref);
		
		
		
		return 1;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
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
	            if ( i == 0 ) System.out.println("WARNING: Indel at the start of the read "+r.getReadName());
	            if ( i == nCigarElems - 1) System.out.println("WARNING: Indel at the end of the read "+r.getReadName());

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
