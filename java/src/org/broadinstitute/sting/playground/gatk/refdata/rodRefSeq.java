package org.broadinstitute.sting.playground.gatk.refdata;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.zip.DataFormatException;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.refdata.BasicReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.xReadLines;

public class rodRefSeq extends BasicReferenceOrderedDatum {

	private GenomeLoc location = null;
	private List<RefSeqRecord> records = null;
	
	public rodRefSeq(String name) {
		super(name);
//		location = new GenomeLoc(0,0,-1);
	}

	public rodRefSeq(String name, GenomeLoc location, List<RefSeqRecord> records) {
		super(name);
		this.location = location;
		this.records = records;
	}

	@Override
	public GenomeLoc getLocation() {
		return location;
	}

	@Override
	public boolean parseLine(Object header, String[] parts) {
		return false; // this rod has its own iterator
	}

	
	public boolean isCoding() {
		if ( records == null ) return false;
		for ( RefSeqRecord r : records) {
			if (location.getStart() >= r.getCodingLocation().getStart() && location.getStart() <= r.getCodingLocation().getStop() ) return true;
		}
		return false;
	}
	

	public boolean isExon() {
		if ( records == null ) return false;
		for ( RefSeqRecord r : records) {
			for ( GenomeLoc e : r.getExons() ) {
				if (location.getStart() >= e.getStart() && location.getStart() <= e.getStop() ) return true;
			}
		}
		return false;
	}

	@Override
	public String repl() {
		throw new StingException("repl() is not implemented yet");
	}

	@Override
	public String toSimpleString() {
		if ( records == null ) return new String(getName()+": <NULL>");
		StringBuilder b = new StringBuilder();
		b.append(getName());
		b.append(":");
		for ( RefSeqRecord r : records ) { 
				b.append(' ');
				b.append(r.getTranscriptId());
		}
		return b.toString();
	}

	@Override
	public String toString() {
		return toSimpleString();
	}

	public static Iterator<rodRefSeq> createIterator(String trackName, File f) throws IOException, DataFormatException {
		return new refSeqIterator(f, trackName);
	}
	
}

class refSeqIterator implements Iterator<rodRefSeq> {
//	private xReadLines reader = null;
	private long curr_position = 0;
	private long max_position = 0;
	private String curr_contig_name = null; // will keep the name of the contig the iterator is currently in
	private List<RefSeqRecord> records; // will keep the list of all transcripts overlapping with the current position
	private PushbackIterator<RefSeqRecord> reader;
	private String name = null;
	
	public refSeqIterator(File f, String trackName) throws IOException {
		reader = new PushbackIterator<RefSeqRecord>( new refSeqRecordIterator(f) );
		records = new LinkedList<RefSeqRecord>();
		name = trackName;
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
		curr_position++;
		if ( curr_position <= max_position ) { 
			// we still have bases covered by at least one currently loaded transcript;
			// we have to purge only subset of transcripts, on which we moved past the end
			Iterator<RefSeqRecord> i = records.iterator();
			while ( i.hasNext() ) {
				RefSeqRecord r = i.next();
				if ( r.getLocation().getStop() < curr_position ) i.remove(); // we moved past the end of transcript r, purge it forever
			}
		} else {
			// ooops, we are past the end of all loaded transcripts - kill them all, 
			// load next transcript and fastforward current position to its start
			records.clear();
			RefSeqRecord r = reader.next(); // if hasNext() previously returned true, we are guaranteed that this call to reader.next() is safe
			records.add( r );  
			curr_contig_name = r.getLocation().getContig();
			curr_position = r.getLocation().getStart();
			max_position = r.getLocation().getStop();
		}
		
		// 'records' only keeps those transcripts now, on which we did not reach the end yet 
		// (we might have reloaded records completely if it was necessary); current position is correctly set.
		// lets check if we walked into new transcripts so we'd need to load them too:
		
		while ( reader.hasNext() ) {
			RefSeqRecord r = reader.peek();
			int ci1 = GenomeLoc.getContigIndex(curr_contig_name);
			int ci2 = GenomeLoc.getContigIndex( r.getLocation().getContig() ); 
			if ( ci1 > ci2 ) throw new StingException("RefSeq track seems to be not contig-ordered");
			if ( ci1 < ci2 ) break; // next transcript is on the next contig, we do not need it yet...
			if ( r.getLocation().getStart() > curr_position ) break; // next transcript is on the same contig but starts after the current position; we are done
			long stop = r.getLocation().getStop();
			if ( stop < curr_position ) throw new StingException("DEBUG: encuntered contig that should have been loaded earlier");
			if ( stop > max_position ) max_position = stop;
			records.add(r);
		}
		
		// 'records' and current position are fully updated. We can now create new rod and return it (NOTE: this iterator will break if the list
		// of pre-loaded records is meddled with by the clients between iterations, so we return them as unmodifiable list)
		rodRefSeq rod = new rodRefSeq(name,new GenomeLoc(curr_contig_name,curr_position, curr_position),Collections.unmodifiableList(records));
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
class refSeqRecordIterator implements Iterator<RefSeqRecord> {
	
	private xReadLines reader;
	
	public refSeqRecordIterator(File f) throws IOException { reader = new xReadLines(f); }

	@Override
	public boolean hasNext() {
		return reader.hasNext();
	}

	@Override
	public RefSeqRecord next()  {
		RefSeqRecord r = new RefSeqRecord();
		r.parseLine( reader.next() );
		return r;
	}
	

	@Override
	public void remove() {
        throw new UnsupportedOperationException();		
	}

}


/** Holds a single transcript annotation: refseq id, gene name, genomic locations of the locus, of the coding region 
 * and of all the exons. 
 */
class RefSeqRecord {
	private String transcript_id;
	private int strand;
	private GenomeLoc transcript_interval;
	private GenomeLoc transcript_coding_interval;
	private List<GenomeLoc> exons;
	private String gene_name;
	private List<Integer> exon_frames;

	public RefSeqRecord() {}
	
	public RefSeqRecord(String line) throws DataFormatException {
		parseLine(line);
	}
	
	public String getTranscriptId() { return transcript_id; }
	public int getStrand() { return strand; }
	public GenomeLoc getLocation() { return transcript_interval; }
	public GenomeLoc getCodingLocation() { return transcript_coding_interval; }
	public String getGeneName() { return gene_name; }
	public int getNumExons() { return exons.size(); }
	public GenomeLoc getExonLocation(int n) {
		if ( n >= exons.size() || n < 0 ) throw new StingException("Index out-of-bounds. Transcript has " + exons.size() +" exons; requested: "+n);
		return exons.get(n);
	}
	public List<GenomeLoc> getExons() { return exons; }
	
	public void parseLine(String line) {
		String[] fields = line.split("\t");
		transcript_id = fields[1];
		if ( fields[3].length()==1 && fields[3].charAt(0)=='+') strand = 1;
		else if ( fields[3].length()==1 && fields[3].charAt(0)=='-') strand = -1;
		else throw new StingException("Expected strand symbol (+/-), found: "+fields[3]);
		
		String contig_name = fields[2];
		transcript_interval = new GenomeLoc(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5]));
		transcript_coding_interval = new GenomeLoc(contig_name, Integer.parseInt(fields[6])+1, Integer.parseInt(fields[7]));
		gene_name = fields[12];
		String[] exon_starts = fields[9].split(",");
		String[] exon_stops = fields[10].split(",");
		String[] eframes = fields[15].split(",");
			
		assert exon_starts.length == exon_stops.length : "Data format error: numbers of exon start and stop positions differ";
		assert exon_starts.length == eframes.length : "Data format error: numbers of exons and exon frameshifts differ";
		
		exons = new ArrayList<GenomeLoc>(exon_starts.length-1);
		exon_frames = new ArrayList<Integer>(eframes.length - 1);
		
		for ( int i = 0 ; i < exon_starts.length - 1 ; i++ ) {
			exons.add(new GenomeLoc(contig_name, Integer.parseInt(exon_starts[i]+1), Integer.parseInt(exon_stops[i]) ) );
			exon_frames.add(Integer.decode(eframes[i]));
		}
	}
	
}