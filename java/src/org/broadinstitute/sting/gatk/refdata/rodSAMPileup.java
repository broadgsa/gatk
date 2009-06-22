package org.broadinstitute.sting.gatk.refdata;


import java.io.IOException;
import java.util.*;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;

import net.sf.picard.reference.ReferenceSequenceFileWalker;

/**
 *  This class wraps Maq/samtools allele calls from pileup format and presents them as a ROD.<br>
 *
 * Example format:<br>
 *   for SNP:<br>
 *      [chr]       [pos] [ref]  [consensus allele(s)]    [consensus confidence]  [snp confidence]   [max mapping qual]  [num reads in the pile] <br>
 *     chrX       466     T                     Y                                   170                                 170                          88                              32 ... (piles of read bases  and quals follow) <br>
 *     <br>
 *   for indel: <br>
 *    [chr]         [pos]     [always *]  [consensus alleles]    [consensus conf.?]        [indel conf.?]     [max mapping qual]     [num reads in the pile]  [indel]   [always *?] [reads haveindel]  [reads may have indel] [other reads?]<br>     
 *     chrX       141444        *           +CA/+CA                               32                                       468            255                                         25                                        +CA     *                         5                                   2                                   12             <br>
 * User: asivache
 * Date: Apr 03, 2009
 * Time: 2:58:33 PM
 * To change this template use File | Settings | File Templates.
 */



public class rodSAMPileup extends BasicReferenceOrderedDatum implements GenotypeList {
	// ----------------------------------------------------------------------
	//
	// Constructors
	//
	// ----------------------------------------------------------------------
	
	private SAMPileupRecord indelGenotype = null;
	private SAMPileupRecord pointGenotype = null;
	
	public rodSAMPileup(final String name) {
		super(name);
	}
	
	@Override
	public GenomeLoc getLocation() {
		if ( pointGenotype != null ) return pointGenotype.getLocation();
		if ( indelGenotype != null ) return indelGenotype.getLocation();
		return null;
	}

	/** Required by ReferenceOrderedDatum interface. This implementation provides its own iterator,
	 * so this method does nothing at all (always returns false).
	 *  
	 */
	@Override
	public boolean parseLine(Object header, String[] parts) throws IOException {
		return false;
	}

	@Override
	public String toString() {
		StringBuilder b = new StringBuilder();
		if ( pointGenotype != null ) {
			b.append(pointGenotype.toString());
			if ( indelGenotype != null ) b.append('\n');
		}
		if ( indelGenotype != null ) b.append(indelGenotype.toString());
		return b.toString();
	}

	@Override
	public Genotype getIndelGenotype() {
		return indelGenotype;
	}

	@Override
	public Genotype getPointGenotype() {
		return pointGenotype;
	}

	@Override
	public boolean hasIndelGenotype() {
		return indelGenotype != null;
	}

	@Override
	public boolean hasPointGenotype() {
		return pointGenotype != null;
	}
	
	/** Iterates over SAM pileup multiline records: each new invocation of next() will 
	 * move to next genomic position in the pileup file, and if the record for that position 
	 * consists of multiple lines (as is the case for indels), they all will be read and processed. 
	 * @author asivache
	 *
	 */
	private static class rodSAMPileupIterator implements Iterator<rodSAMPileup> {
		//private xReadLines parser = null;
		private String rodName = null;
		private PushbackIterator<SAMPileupRecord> parser = null;
		private long line = 0;
		
		rodSAMPileupIterator(String name, java.io.File f) {
			parser = new PushbackIterator<SAMPileupRecord>(SAMPileupRecord.createIterator(name,f));
			rodName = name;
		}
		
		@Override
		public boolean hasNext() {
			return parser.hasNext();
		}

		@Override
		public rodSAMPileup next() {

			rodSAMPileup result = new rodSAMPileup(rodName);
			
			SAMPileupRecord r = parser.next();
			
	//		if ( (++line)%10000 == 0 ) System.out.printf("%s: record is %d (%s)%n", rodName, line,r.getLocation().toString());

			if ( r.isPointGenotype() ) result.pointGenotype = r;
			else result.indelGenotype = r;
			
			if ( parser.hasNext() ) {

				SAMPileupRecord r2 = parser.next();
				int cmp = r.getLocation().compareTo(r2.getLocation());
				switch ( cmp ) {
				case -1 : // next record is at greater position
					parser.pushback(r2); // return next record back to the stream
					break;
				case 0: // next record is at the same position; in this case r MUST be point and r2 MUST be indel
					if ( result.indelGenotype != null ) { // oops, we've seen indel already
						if ( r2.isPointGenotype() ) throw new StingException("SAM pileup file might be corrupted: point genotype follows indel genotype line at "+r.getLocation() );
						else throw new StingException("SAM pileup file might be corrupted: two indel genotype lines found at same position "+r.getLocation() );
					} else {
						// ok, the first genotype we've seen was a point one
						if ( r2.isPointGenotype() ) throw new StingException("SAM pileup file might be corrupted: two point genotype lines found at same position "+r.getLocation() );
						// and the second is an indel, wheww...
						result.indelGenotype = r2;
					}
					
					// the following is solely for the purposes of strict validation of pileup file: since we've already read two lines at the
					// the same genomic location (point and indel variants), there can not be any more:
					if ( parser.hasNext() ) {
						r2 = parser.next();
						cmp = r.getLocation().compareTo(r2.getLocation() ) ;
						if ( cmp == 0 ) throw new StingException("SAM pileup file is corrupted: more than two lines found at position "+r.getLocation());
						if ( cmp > 0 ) throw new StingException("SAM pileup file or reference dictionary is corrupted: lesser position "+r2.getLocation()+" is encountered right after position " 
								+ r.getLocation() ); 
						parser.pushback(r2);
					}
					
					break;
				case 1: // next record is at lower position?? come'on!
					throw new StingException("SAM pileup file or reference dictionary is corrupted: lesser position "+r2.getLocation()+" is encountered right after position " + r.getLocation() );
				default:
					throw new StingException("INTERNAL ERROR: this point should never be reached");
				}
			}
			return result ;
		}


		@Override
		public void remove() {
			throw new UnsupportedOperationException("'remove' operation is not supported for file-backed SAM pileups");
		}
		
	}
	
	public static Iterator<rodSAMPileup> createIterator(String name, java.io.File file) {
		return new rodSAMPileup.rodSAMPileupIterator(name,file);
	}

	public static void main(String argv[]) {
		String testFile = "/humgen/gsa-scr1/asivache/TCGA/Ovarian/C2K/0805/normal.pileup";
//		String testFile = "/humgen/gsa-scr1/asivache/trios/CEU/NA12891.12892.12878/mother.chr1.pileup.indel";
		
		Iterator<rodSAMPileup> it = createIterator("test-normal", new java.io.File(testFile));
		
		ReferenceSequenceFileWalker reference = new ReferenceSequenceFileWalker(
              new java.io.File(  "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
//              new java.io.File(  "/humgen/gsa-scr1/asivache/trios/CEU/NA12891.12892.12878/human_b36_both.fasta")
                     );

		if ( reference.getSequenceDictionary() == null ) {
			System.out.println("No reference sequence dictionary found. Abort.");
			System.exit(1);
		}

		GenomeLocParser.setupRefContigOrdering(reference.getSequenceDictionary());
		
		int counter = 0;
		
		while ( it.hasNext() && counter < 430 ) {
			rodSAMPileup p = it.next();
			System.out.println(p.toString());
/*			
			System.out.print(p.getLocation().toString());
			System.out.print('\t');
			
			if ( p.isIndel() && p.isSNP() ) { System.out.print("Indel+SNP"); }
			else {
				if ( p.isSNP() ) { System.out.print("SNP"); }
				else { 
					if ( p.isIndel() ) { System.out.print("Indel"); }
					else { System.out.print("REF"); }
				}
			}
			
			System.out.print('\t');
			System.out.print(p.getFWDAlleles().get(0)+"/"+p.getFWDAlleles().get(1));
			System.out.print('\t');
			System.out.println(p.getConsensusConfidence()+"\t"+p.getVariantConfidence());
			*/
			counter++;
		}
	}


}


