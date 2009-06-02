package org.broadinstitute.sting.gatk.refdata;

import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

/** Holds a single transcript annotation: refseq id, gene name, genomic locations of the locus, of the coding region 
 * and of all the exons. 
 */
public class Transcript {
		private String transcript_id;
		private int strand;
		private GenomeLoc transcript_interval;
		private GenomeLoc transcript_coding_interval;
		private List<GenomeLoc> exons;
		private String gene_name;
		private List<Integer> exon_frames;

		public Transcript() {}
		
		public Transcript(String line) {
			parseLine(line);
		}
		
		/** Returns id of the transcript (RefSeq NM_* id) */
		public String getTranscriptId() { return transcript_id; }
		/** Returns coding strand of the transcript, 1 or -1 for positive or negative strand, respectively */
		public int getStrand() { return strand; }
		/** Returns transcript's full genomic interval (includes all exons with UTRs) */  
		public GenomeLoc getLocation() { return transcript_interval; }
		/** Returns genomic interval of the coding sequence (does not include UTRs, but still includes introns, since it's a single interval on the DNA) */  
		public GenomeLoc getCodingLocation() { return transcript_coding_interval; }
		/** Name of the gene this transcript corresponds to (NOT gene id such as Entrez etc) */
		public String getGeneName() { return gene_name; }
		/** Number of exons in this transcript */
		public int getNumExons() { return exons.size(); }
		/** Genomic location of the n-th exon; throws an exception if n is out of bounds */
		public GenomeLoc getExonLocation(int n) {
			if ( n >= exons.size() || n < 0 ) throw new StingException("Index out-of-bounds. Transcript has " + exons.size() +" exons; requested: "+n);
			return exons.get(n);
		}
		/** Returns the list of all exons in this transcript, as genomic intervals */
		public List<GenomeLoc> getExons() { return exons; }
		
		/** Returns true if the specified interval 'that' overlaps with the full genomic interval of this transcript */ 
	    public boolean overlapsP (GenomeLoc that) {
	    	return transcript_interval.overlapsP(that);
	    }
		
		/** Returns true if the specified interval 'that' overlaps with the coding genomic interval of this transcript. 
		 * NOTE: since "coding interval" is still a single genomic interval, it will not contain UTRs of the outermost exons,
		 * but it will still contain introns and/or exons internal to this genomic locus that are not spliced into this transcript. 
		 * @see overlapsExonP()
		 */ 
	    public boolean overlapsCodingP (GenomeLoc that) {
	    	return transcript_coding_interval.overlapsP(that);
	    }

		/** Returns true if the specified interval 'that' overlaps with any of the exons actually spliced into this transcript */ 
	    public boolean overlapsExonP (GenomeLoc that) {
	    	for ( GenomeLoc e : exons ) {
	    		if ( e.overlapsP(that) ) return true;
	    	}
	    	return false;
	    }

	    /** Fills this object from a text line in RefSeq (UCSC) text dump file */
		public void parseLine(String line) {
			String[] fields = line.split("\t");
			transcript_id = fields[1];
			if ( fields[3].length()==1 && fields[3].charAt(0)=='+') strand = 1;
			else if ( fields[3].length()==1 && fields[3].charAt(0)=='-') strand = -1;
			else throw new StingException("Expected strand symbol (+/-), found: "+fields[3]);
			
			String contig_name = fields[2];
			transcript_interval = GenomeLoc.parseGenomeLoc(contig_name, Integer.parseInt(fields[4])+1, Integer.parseInt(fields[5]));
			transcript_coding_interval = GenomeLoc.parseGenomeLoc(contig_name, Integer.parseInt(fields[6])+1, Integer.parseInt(fields[7]));
			gene_name = fields[12];
			String[] exon_starts = fields[9].split(",");
			String[] exon_stops = fields[10].split(",");
			String[] eframes = fields[15].split(",");
				
			assert exon_starts.length == exon_stops.length : "Data format error: numbers of exon start and stop positions differ";
			assert exon_starts.length == eframes.length : "Data format error: numbers of exons and exon frameshifts differ";
			
			exons = new ArrayList<GenomeLoc>(exon_starts.length);
			exon_frames = new ArrayList<Integer>(eframes.length);
			
			for ( int i = 0 ; i < exon_starts.length  ; i++ ) {
				exons.add(GenomeLoc.parseGenomeLoc(contig_name, Integer.parseInt(exon_starts[i])+1, Integer.parseInt(exon_stops[i]) ) );
				exon_frames.add(Integer.decode(eframes[i]));
			}
		}
		
}

