package org.broadinstitute.sting.gatk.refdata;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.refdata.BasicReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.Transcript;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.xReadLines;

public class rodRefSeq extends BasicReferenceOrderedDatum {

	private GenomeLoc location = null;
	private List<Transcript> records = null;
	
	public rodRefSeq(String name) {
		super(name);
//		location = GenomeLoc.parseGenomeLoc(0,0,-1);
	}

	/** Despite this constructor is public, it is meant primarily for the internal use; RefSeq iterator will
	 * call it to populate the ROD at the given genomic location with the data (transcripts) it is holding
	 * @param name
	 * @param location
	 * @param records
	 */
	public rodRefSeq(String name, GenomeLoc location, List<Transcript> records) {
		super(name);
		this.location = location;
		this.records = records;
	}

	@Override
	public GenomeLoc getLocation() {
		return location;
	}

	/** Required by ReferenceOrderedDatum interface; this method does nothing (always returns false),
	 * since this rod provides its own iterator for reading underlying data files.
	 */
	@Override
	public boolean parseLine(Object header, String[] parts) {
		return false; // this rod has its own iterator
	}

	/** Returns true if the current position this ROD is associated with is within the coding interval for at least
	 * one of the annotated transcripts. NOTE: "coding" interval is defined as a single genomic interval, so it
	 * does not include the UTRs of the outermost exons, but it includes introns between exons spliced into a 
	 * transcript, or internal exons that are not spliced into a given transcript. To check that a position is
	 * indeed within an exon but not in UTR, use isExon() && isCoding(). @see isExon().
	 * @return
	 */
	public boolean isCoding() {
		if ( records == null ) return false;
		for ( Transcript t : records) {
			if ( t.overlapsCodingP(location) ) return true;
		}
		return false;
	}
	

	/** Returns true if the current position this ROD is associated with is within an exon for at least
	 * one of the annotated transcripts. NOTE: position can be still within a UTR, see @isCoding()
	 * @return
	 */
	public boolean isExon() {
		if ( records == null ) return false; 
		for ( Transcript t : records) {
			if ( t.overlapsExonP(location) ) return true;
		}
		return false;
	}

	/** Returns all annotated transcripts overlapping with the current position as an immutable list.
	 * "Overlap" is defined as position being within genomic interval corresponding to the whole 
	 * transcript start/transcript stop coordinates, thus the position can be still in a UTR, in intron, or
	 * in an internal exon that is actually not spliced into the specific transcript. Use isExon(), isCoding() of
	 * this ROD class, or query individual Transcript objects returned by this method to get more details. 
	 * @return
	 */
	public List<Transcript> getTranscripts() { return Collections.unmodifiableList(records) ; }
	
	@Override
	public String repl() {
		throw new StingException("repl() is not implemented yet");
	}

	/** Will print the genomic location of this rod, followed by (space separated) ids of all the
	 * annotated transcripts overlapping with this position. 
	 */
	@Override
	public String toSimpleString() {
		if ( records == null ) return new String(getName()+": <NULL>");
		StringBuilder b = new StringBuilder();
		b.append(getName());
		b.append(":");
		for ( Transcript t : records ) { 
				b.append(' ');
				b.append(t.getTranscriptId());
		}
		return b.toString();
	}

	@Override
	public String toString() {
		return toSimpleString();
	}

	public static Iterator<rodRefSeq> createIterator(String trackName, File f) throws IOException {
//		System.out.println("REFSEQ ITERATOR CREATED");
		return new refSeqIterator(trackName,f);
	}
	
}

class refSeqIterator implements Iterator<rodRefSeq> {
//	private xReadLines reader = null;
	private long curr_position = 0;
	private long max_position = 0;
	private String curr_contig_name = null; // will keep the name of the contig the iterator is currently in
	private List<Transcript> records; // will keep the list of all transcripts overlapping with the current position
	private PushbackIterator<Transcript> reader;
	private String name = null;
//	private long z = 0; // counter, for debugging only
//	private long t = 0; // for debug timer only
	
	public refSeqIterator(String trackName, File f) throws IOException {
		reader = new PushbackIterator<Transcript>( new refSeqRecordIterator(f) );
		records = new LinkedList<Transcript>();
		name = trackName;
//		System.out.println("REFSEQ ITERATOR CONSTRUCTOR");
	}

	@Override
	public boolean hasNext() {
		// if we did not walk to the very end of currently loaded transcripts, then we 
		// definitely have next
		if ( curr_position < max_position ) return true;
		
		// we are past currently loaded stuff; we have next if there are more lines to load:
		return reader.hasNext();
	}

	@Override
	public rodRefSeq next() {
//		if ( z == 0 ) t = System.currentTimeMillis();
		curr_position++;
		if ( curr_position <= max_position ) { 
			// we still have bases covered by at least one currently loaded transcript;
			// we have to purge only subset of transcripts, on which we moved past the end
			Iterator<Transcript> i = records.iterator();
			while ( i.hasNext() ) {
				Transcript t = i.next();
				if ( t.getLocation().getStop() < curr_position ) {
					i.remove(); // we moved past the end of transcript r, purge it forever
				}
			}
		} else {
			// ooops, we are past the end of all loaded transcripts - kill them all, 
			// load next transcript and fastforward current position to its start
			records.clear();
			Transcript t = reader.next(); // if hasNext() previously returned true, we are guaranteed that this call to reader.next() is safe
			records.add( t );  
			curr_contig_name = t.getLocation().getContig();
			curr_position = t.getLocation().getStart();
			max_position = t.getLocation().getStop();
		}
		
		// 'records' only keeps those transcripts now, on which we did not reach the end yet 
		// (we might have reloaded records completely if it was necessary); current position is correctly set.
		// lets check if we walked into additional new transcripts so we'd need to load them too:
		
		while ( reader.hasNext() ) {
			Transcript t = reader.peek();
			int ci1 = GenomeLoc.getContigIndex(curr_contig_name);
			int ci2 = GenomeLoc.getContigIndex( t.getLocation().getContig() ); 
			if ( ci1 > ci2 ) throw new StingException("RefSeq track seems to be not contig-ordered");
			if ( ci1 < ci2 ) break; // next transcript is on the next contig, we do not need it yet...
			if ( t.getLocation().getStart() > curr_position ) break; // next transcript is on the same contig but starts after the current position; we are done
			t = reader.next(); // we do need next record, time to load it for real
			long stop = t.getLocation().getStop();
			if ( stop < curr_position ) throw new StingException("DEBUG: encountered contig that should have been loaded earlier");
			if ( stop > max_position ) max_position = stop;
			records.add(t);
		}
		
		// 'records' and current position are fully updated. We can now create new rod and return it (NOTE: this iterator will break if the list
		// of pre-loaded records is meddled with by the clients between iterations, so we return them as unmodifiable list)
		rodRefSeq rod = new rodRefSeq(name,GenomeLoc.parseGenomeLoc(curr_contig_name,curr_position, curr_position),Collections.unmodifiableList(records));
//		if ( (++z) % 1000000 == 0 ) {
//			System.out.println(rod.getLocation()+": holding "+records.size()+ "; time per 1M ref positions: "+((double)(System.currentTimeMillis()-t)/1000.0)+" s");
//			z = 0;
//		}
		return rod;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();		
	}
	
	

}

/** Low-level iterator for reading refseq annotation file record by record (i.e. line by line). Returns 
 * pre-processed input lines as RefSeqRecord objects.
 */
class refSeqRecordIterator implements Iterator<Transcript> {
	
	private xReadLines reader;
	
	public refSeqRecordIterator(File f) throws IOException { reader = new xReadLines(f); }

	@Override
	public boolean hasNext() {
		return reader.hasNext();
	}

	@Override
	public Transcript next()  {
		Transcript t = new Transcript();
		t.parseLine( reader.next() );
		return t;
	}
	

	@Override
	public void remove() {
        throw new UnsupportedOperationException();		
	}

}


