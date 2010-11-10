/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.playground.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.*;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.collections.Pair;

public class GenomicMap implements Iterable<Map.Entry<String, Collection<GenomeLoc> > >{
	
	private Map<String,Collection<GenomeLoc> > map;
	
	/** Creates new empty genomic map preallocated to handle <code>initialContig</code> contigs.
	 * 
	 * @param initialContigs
	 */
	public GenomicMap(int initialContigs) {
		map = new HashMap<String,Collection<GenomeLoc> >(initialContigs);
	}
	
	/** Creates new empty genomic map */
	public GenomicMap() {
		this(1000);
	}

	/** Adds custom contig to the map, as a collection of intervals on the master reference.
	 * 
	 * @param name name of the custom contig; can not be null
	 * @param c mapping of the custom contig sequence onto intervals on the master reference
	 */
	public void addCustomContig(String name, Collection<GenomeLoc> c) {
		if ( name == null ) throw new ReviewedStingException("Custom contig name can not be null");
		if ( map.containsKey(name)) throw new ReviewedStingException("Custom contig "+name+" already exists");
		map.put(name, c);
	}
	
	/** Returns mapping of the specified custom contig onto the custom reference. \
	 * If no such contig exists, returns null.
	 * 
	 * @param name
	 * @return
	 */
	public Collection<GenomeLoc> getContigMapping(String name) { return map.get(name); }
	
	/** Read genomic map from specified Arachne multimap file. Format: 
	 * contig_id start stop contig_id start stop ... # name ...
	 * where start, stop are 0 based, closed intervals  
	 * @param f
	 */
	public void readArachne(SAMSequenceDictionary sequenceDictionary,GenomeLocParser genomeLocParser,File f) {
		
		try {
			BufferedReader reader = new BufferedReader( new FileReader(f) );
			
			String line = null;
			while( ( line = reader.readLine() ) != null ) {
				String[] halves = line.split("#",2);
				if ( halves.length < 2 ) 
					throw new UserException.MalformedFile(f, "Line: "+line+"\nin map file "+f+"\n does not contain contig name");
				
				int p1 = 0;
				for ( ;  p1 < halves[1].length() && Character.isWhitespace(halves[1].charAt(p1) ); p1++ ); 
				// p1 is now index of first non-space
				int p2 = p1;
				for ( ; p2 < halves[1].length() && ! Character.isWhitespace(halves[1].charAt(p2) ); p2++ );
				// p2 is index of first whitespace after first word
				
				if ( p1 == p2 ) 
					throw new UserException.MalformedFile(f, "Line: "+line+"\n in map file "+f+"\nNo contig name found after '#'");
				
				String name = halves[1].substring(p1, p2);
								
				String[] coord_parts = halves[0].split("\\s");
				if ( coord_parts.length % 3 != 0 ) 
					throw new UserException.MalformedFile(f, "Line: "+line+"\n in map file "+f+"\nNumber of coordinate fields is not a multiple of 3");
				
				List<GenomeLoc> segments = new ArrayList<GenomeLoc>( coord_parts.length / 3 );
				
				for ( int i = 0 ; i < coord_parts.length ; i += 3 ) {
					// Arachne map file contains 0-based, closed intervals, hence +1 below.
					int index = Integer.parseInt(coord_parts[i]);
                    String contig = sequenceDictionary.getSequence(index).getSequenceName();
					int start = Integer.parseInt(coord_parts[i+1]);
					int stop = Integer.parseInt(coord_parts[i+2]);
					segments.add(genomeLocParser.createGenomeLoc(contig, start+1, stop+1));
				}
				
				addCustomContig(name, segments);
				
			}
			reader.close();
		} catch ( FileNotFoundException e) {
			throw new UserException.CouldNotReadInputFile(f, e);
		} catch (IOException e) {
			throw new UserException.CouldNotReadInputFile(f, e);
		}
	}

	/** Read genomic map from specified file in "new" format. Format: 
	 * name chr:start-stop,chr:start-stop,...,chr:start-stop
	 * where start, stop are 1 based, closed intervals  
	 * @param f
	 */
	public void read(GenomeLocParser genomeLocParser,File f) {
		
		try {
			BufferedReader reader = new BufferedReader( new FileReader(f) );

			String line = null;
			while( ( line = reader.readLine() ) != null ) {
				int p1 = 0;
				while ( p1 < line.length() && Character.isWhitespace(line.charAt(p1))) p1++;
				int p2 = p1;
				while ( p2 < line.length() && ! Character.isWhitespace(line.charAt(p2))) p2++;
				if ( p1 == p2 ) continue; // empty line

				String name = line.substring(p1, p2);
				
				List<GenomeLoc> segments = new ArrayList<GenomeLoc>( 5 );

				p1 = p2+1; // set p1 after first whitespace after the name
				while ( p1 < line.length() && Character.isWhitespace(line.charAt(p1))) p1++; // skip whitespaces
				p2 = p1;
				while ( p2 < line.length() && line.charAt(p2) != ',') p2++; // next comma or end-of-line

				while ( p2 != p1 ) {
					GenomeLoc newSegment = genomeLocParser.parseGenomeLoc(line.substring(p1, p2));
					if ( segments.size() > 0 &&
							segments.get(segments.size()-1).getStop()+1 == newSegment.getStart() &&
							segments.get(segments.size()-1).getContigIndex() == newSegment.getContigIndex())
						System.out.println("WARNING: strictly adjacent segments found in custom contig "+name);

					segments.add(newSegment);
				
					p1 = p2+1; // set p1 after the comma
					while ( p1 < line.length() && Character.isWhitespace(line.charAt(p1))) p1++; // skip whitespaces
					p2 = p1;
					while ( p2 < line.length() && line.charAt(p2) != ',') p2++; // next comma or end-of-line
				}
				if ( segments.size() == 0 ) throw new ReviewedStingException("Line "+line+" has no intervals specified");
				addCustomContig(name, segments);
			}
			reader.close();
		} catch ( FileNotFoundException e) {
			throw new UserException.CouldNotReadInputFile(f, e);
		} catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(f, e);
		}
	}

	public void write(File f) {
		try {
			BufferedWriter writer = new BufferedWriter( new FileWriter( f ));
			for ( String name : nameSet() ) {
				writer.append(name+" ");
				Iterator<GenomeLoc> iter = getContigMapping(name).iterator();
				if ( iter.hasNext() ) writer.append(iter.next().toString());
				while (iter.hasNext()) {
					writer.append(',');
					writer.append(iter.next().toString());
				}
				writer.append('\n');
			}
			writer.close();
		} catch (IOException e) {
			throw new UserException.CouldNotCreateOutputFile(f, e);
		}
	}
	
	/** Remaps a record (read) aligned to a custom contig back onto the master reference. 
	 * If the map does not have mapping information for
	 * the contig, an exception will be thrown. This method changes read's reference name, start position and 
	 * cigar, as well as the read's file header (must be provided). 
	 * 
	 * Some aligners (e.g. bwa) can return "alignments" spanning across contig boundaries. The last argument of this
	 * method controls the behavior in this case: if it is set to true, such alignments are ignored upon detection,
	 * and the method returns null. Otherwise, strict validation mode is used: if aligned read extends beyond the
	 * contig boundary, an exception is thrown.
	 *  
	 * @param r read, alignment information (contig, start position, cigar) will be modified by this method
	 * @param h SAM file header for the master reference the alignment is being mapped onto; will be substituted for the read's header.
	 * @return same read instance that was passed to this method, remapped
	 */
	public SAMRecord remapToMasterReference(SAMRecord r, SAMFileHeader h, boolean discardCrossContig) {
		if ( AlignmentUtils.isReadUnmapped(r) ) {
            // set to NO_... just in case: in principle, SAM format spec allows unmapped reads (with 'unmapped'
            // flag raised) to have reference contig and start position set to arbitrary values for sorting
            // purposes; after remapping, these values would make no sense or even cause a crash when reading
            // remapped bam
            r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            return r; // nothing to do if read is unmapped
        }
		
		int customStart = r.getAlignmentStart();

        // get mapping from read's contig onto a "global" contig (as a list of intervals on the latter):
		Collection<GenomeLoc> segments = getContigMapping(r.getReferenceName());
		if ( segments == null ) throw new UserException.MalformedBam(r, "Can not remap a record: unknown custom contig name "+r.getReferenceName());

        // scroll the list of intervals until we find the interval that the alignment start falls into:
		Pair<? extends Iterator<GenomeLoc>, Integer> p = seekForward(segments,customStart);
		
		Iterator<GenomeLoc> iter = p.first;
		
		GenomeLoc gl = iter.next(); // initialization: get interval that contains customStart

        // p.second is 1-based position of alignment start relative to interval gl;
        // hence refPos is 1-based position of the alignment start on the master ref (since gl.getStart() is 1-based too)
		int refPos = (int)(p.second+gl.getStart()-1); 
		
		String oldRefName = r.getReferenceName();
		int oldStart = r.getAlignmentStart();
		int oldEnd = r.getAlignmentEnd();

		r.setAlignmentStart(refPos);		
		
		r.setHeader(h); // have to substitute here, or setReferenceIndex will not work correctly below
		
		r.setReferenceIndex(gl.getContigIndex());
		
		Cigar oldCigar = r.getCigar();
		Cigar newCigar = new Cigar();
		int N = oldCigar.numCigarElements() ;
		
		long currStop = gl.getStop();// end of the current segment of the custom contig on the master reference, 1-based inclusive
        int delayedGap = 0 ; // length of the 'N' gap between the segments (intervals) on the master ref, to be added only if followed by another cigar element		
		for ( int k = 0; k < N ; k++ ) {
			CigarElement ce = oldCigar.getCigarElement(k);
			int len = ce.getLength();      // length of the cigar element
			switch( ce.getOperator() ) {
			case S: // soft clip
			case H: // or hard clip - these are not included in getAlignmentStart, so pass them through
				if ( k != 0 && k != N-1 ) // paranoid
					throw new ReviewedStingException("Don't know what to do with S or N cigar element that is not at the either end of the cigar. Cigar: "+
							r.getCigarString());
			case I: // insertions are passed through as well
				newCigar.add(new CigarElement(len,ce.getOperator()));
				break;
			case D:
            case M:
///////////
                if ( delayedGap > 0 ) {
                    // we get here if previous M or D element ended exactly at the interval boundary; we need
                    // to add the stretch of N's only if that element turned out to be not the last one, so we do it now
                    newCigar.add(new CigarElement(delayedGap, CigarOperator.N));
                    delayedGap = 0;
                }
                while ( refPos + len - 1 > currStop ) { // current D or M cigar element extends beyond the end of current segment

                    // we have that many bases in the current cigar element till the end of the current segment:
                    int currLength = (int)(currStop-refPos+1);


                    // curr length can be exactly 0 if previous element ended exactly at the segment boundary:
                    // after that element was processed, refPos was set to currStop+1, so in this special case we need
                    // *first* to switch to next segment, *then* start adding bases from the current element.
                    if ( currLength > 0 ) {
                        newCigar.add(new CigarElement( currLength,ce.getOperator()) ); // record deletion/match till the end of the current segment
                        len -= currLength; // we still have 'len' bases remaining in the current cigar element
                    }

                    // NOTE: since we entered the loop, we were guaranteed that len > currLength, so now len > 0

                    // check if we have next segment to extend remaining matching bases to; if we don't, something's awfully wrong:
                    if ( ! iter.hasNext() ) {
                        String message = "Record "+r.getReadName()+" extends beyond its custom contig."+
                                        "\nRead aligns to: "+oldRefName+":"+oldStart+"-"+oldEnd+"; cigar="+
                                        r.getCigarString()+"; contig length="+contigLength(segments);
                        if ( discardCrossContig ) {
                //			System.out.println("WARNING: ALIGNMENT DISCARDED: "+message);
                            return null;
                        } else throw new UserException.MalformedBam(r, message);
                    }

                    gl = iter.next(); // advance to next segment

                    refPos = (int)gl.getStart(); // we jump to the start of next segment on the master ref

                    if ( gl.getContigIndex() != r.getReferenceIndex() )
                        throw new UserException.MalformedBam(r, "Contig "+oldRefName+
                        " has segments on different master contigs: currently unsupported");

                    if ( refPos < currStop + 1 )
                        throw new UserException.MalformedBam(r, "Contig "+oldRefName+
                        " has segments that are out of order or strictly adjacent: currently unsupported");
                    if ( len > 0 && refPos > currStop + 1 ) {
                        // add "panning" N's w/respect to the master ref over the region between adjacent segments
                        // (and do not add anything if segments are strictly adjacent, i.e. refPos == currStop+1):
                        newCigar.add(new CigarElement((int)(refPos-currStop-1),CigarOperator.N));
                    } else {
                        // we jumped onto the next interval, but the current cigar element ended exactly
                        // at the end of the previous interval. We will need to end later, only if more M/D elements follow:
                        delayedGap = (int)(refPos-currStop-1);
                    }
                    currStop = gl.getStop();
                    // now we can continue with recording remaining matching bases over the current segment
                }
                // we get here when remaining matching bases fit completely inside the current segment:
                if ( len > 0 ) newCigar.add(new CigarElement(len,ce.getOperator()));
                refPos+=len;

                break;
////////////
			}
		}
		
		r.setCigar(newCigar);
		
		return r;
	}
	
	public int size() { return map.size(); }
	
	public Iterator<Map.Entry<String, Collection<GenomeLoc> > > iterator() { return map.entrySet().iterator(); }
	public Iterator<String> nameIterator() { return map.keySet().iterator(); }
	public Set<String> nameSet() { return map.keySet(); }
	
	/** Returns an iterator into the specified collection of segments that points right before the segment that contains
	 * specified position, and the offset of the position inside that segment. This helper method assumes that
	 * there is a "custom" contig built of intervals on the "master" reference; the first argument specifies
	 * the mapping (i.e. an ordered collection of master reference intervals the custom contig is built of), and the second argument
	 * is the 1-based position on that custom contig. Returned iterator is advanced towards the interval (element of the passed
	 * collection) that contains the specified position, namely a call to next() on the returned iterator will return that interval.
	 * Returned integer offset is the 1-based offset of the base at position <code>position</code> on the custom contig with respect
	 * to the start of the interval that base. If position is outside of the custom contig, runtime StingException will be thrown. 
	 * @param segments mapping of the custom contig onto the master reference
	 * @param position 1-based position on the custom contig
	 * @return
	 */
	private Pair<PushbackIterator<GenomeLoc>,Integer> seekForward(Collection<GenomeLoc> segments,int position) {
		
		if ( position < 1 ) throw new ReviewedStingException("Position "+position + " is outside of custom contig boundaries");
		
		PushbackIterator<GenomeLoc> iter = new PushbackIterator<GenomeLoc>(segments.iterator());
		
		while ( iter.hasNext() ) {
			GenomeLoc current = iter.next();
			long length = current.getStop() - current.getStart() + 1; // length of current segment
			if ( position <= length ) { // position is on the current segment
				iter.pushback(current);
				return new Pair<PushbackIterator<GenomeLoc>, Integer >( iter,position);
			}
			// no, position is beyond the current segment; subtract the length of current segment and step to next one
			position -= length;
		}
		// if we get here, position is to the right of the last segment; not good.
		throw new ReviewedStingException("Position "+position + " is outside of custom contig boundaries");
	}

	private long contigLength(Collection<GenomeLoc> segments) {
		long l = 0;
		for ( GenomeLoc g : segments ) l += (g.getStop() - g.getStart() + 1 );
		return l;
	}
	
	public static void main(String argv[]) {
		
//		SAMFileReader reader = new SAMFileReader(new java.io.File("/humgen/gsa-scr1/asivache/TCGA/Ovarian/C2K/0904/normal.bam"));
        SAMFileReader reader = new SAMFileReader(new java.io.File("X:/asivache/cDNA/new_pipeline/30BV1/test.1.sam"));

		
		SAMRecord r = new SAMRecord(reader.getFileHeader());
        GenomeLocParser genomeLocParser = new GenomeLocParser(reader.getFileHeader().getSequenceDictionary());

        r.setReferenceName("ENST00000378466");
        r.setAlignmentStart(1235);
        r.setCigarString("24M1D27M");

//		List<GenomeLoc> s = new ArrayList<GenomeLoc>();
//		s.add( GenomeLocParser.createGenomeLoc("chr1", 100, 199));	
//		s.add( GenomeLocParser.createGenomeLoc("chr1", 300, 499));	
//		s.add( GenomeLocParser.createGenomeLoc("chr1", 600, 799));	

		GenomicMap m = new GenomicMap(5);
		
//		m.readArachne(genomeLocParser,new File("/humgen/gsa-scr1/asivache/cDNA/Ensembl48.transcriptome.map"));
//        m.write(new File("/humgen/gsa-scr1/asivache/cDNA/new_pipeline/Ensembl48.new.transcriptome.map"));
        m.read(genomeLocParser,new File("W:/berger/cDNA_BAM/refs/Ensembl52.plus.Genome.map"));

        m.remapToMasterReference(r,reader.getFileHeader(),true);

//		if ( m.getContigMapping("ENST00000302418") == null ) System.out.println("ERROR! CONTIG IS MISSING!");

		int cnt = 0;
		
		System.out.println(m.size() + " contigs loaded");
        System.out.println("new alignment: "+r.format())  ;
/*
 		for ( String name : m.nameSet() ) {
 
			System.out.print(name);
			System.out.print(": ");
			for ( GenomeLoc g : m.getContigMapping(name)) {
				System.out.print(g.toString()+",  ");
			}
			System.out.println();
			cnt ++;
			if ( cnt > 10 ) break;
		}
*/		
//		m.addCustomContig("My", s);
/*		
		r.setReferenceName("My");
		r.setAlignmentStart(3);
		r.setCigarString("5S97M5D197M5H");
		
		m.remapToMasterReference(r);
		System.out.println(r.getReferenceName()+":"+r.getAlignmentStart()+" "+r.getCigarString());
*/		
		reader.close();
			
	}
	
}
