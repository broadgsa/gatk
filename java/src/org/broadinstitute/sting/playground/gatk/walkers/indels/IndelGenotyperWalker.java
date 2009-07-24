package org.broadinstitute.sting.playground.gatk.walkers.indels;

import net.sf.samtools.*;

import org.broadinstitute.sting.gatk.refdata.RODIterator;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.Transcript;
import org.broadinstitute.sting.gatk.refdata.rodRefSeq;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.playground.utils.CircularArray;
import org.broadinstitute.sting.utils.GenomeLoc;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;


@ReadFilters({Platform454Filter.class, ZeroMappingQualityReadFilter.class})
public class IndelGenotyperWalker extends ReadWalker<Integer,Integer> {
    @Argument(fullName="outputFile", shortName="O", doc="output file name (defaults to BED format)", required=true)
    java.io.File bed_file;
    @Argument(fullName="1kg_format", shortName="1kg", doc="output in 1000 genomes format", required=false)
    boolean FORMAT_1KG;
	@Argument(fullName="somatic", shortName="somatic",
			doc="Perform somatic calls; two input alignment files must be specified", required=false)
	boolean call_somatic = false;
	@Argument(fullName="verbose", shortName="verbose", 
			doc="Tell us what you are calling now (printed to stdout)", required=false)
	boolean verbose = false;
	@Argument(fullName="minCoverage", shortName="minCoverage", 
			doc="must have minCoverage or more reads to call indel; with --somatic this value is applied to tumor sample", required=false)
	int minCoverage = 6;
	@Argument(fullName="minNormalCoverage", shortName="minNormalCoverage", 
			doc="used only with --somatic;  normal sample must have at least minNormalCoverage or more reads to call germline/somatic indel", required=false)
	int minNormalCoverage = 4;
	@Argument(fullName="minFraction", shortName="minFraction", 
			doc="minimum fraction of reads with CONSENSUS indel at a site, out of all reads covering the site, required for a consensus call"+
			" (fraction of non-consensus indels at the site is not considered here, see minConsensusFraction)", required=false)
	double minFraction = 0.3;
	@Argument(fullName="minConsensusFraction", shortName="minConsensusFraction", 
			doc="Minimum fraction of CONSENSUS indel observations at a site wrt all indel observations at the site required to make the call", required=false)
	double minConsensusFraction = 0.7;
	@Argument(fullName="refseq", shortName="refseq", 
			doc="Name of RefSeq transcript annotation file. If specified, indels will be annotated as GENOMIC/UTR/INTRON/CODING", required=false)
	String RefseqFileName = null;
	 
	private static int WINDOW_SIZE = 200;
	private RunningCoverage coverage;
	private RunningCoverage normal_coverage; // when performing somatic calls, we will be using this one for normal, and 'coverage' for tumor
	private int currentContigIndex = -1;
	private int currentPosition = -1; // position of the last read we've seen on the current contig
	private String refName = null;
	private java.io.Writer output = null;
	private GenomeLoc location = null;
	
	private RODIterator<rodRefSeq> refseqIterator=null;

	private Set<String> normalReadGroups;
	private Set<String> tumorReadGroups ;
	
	private int MISMATCH_WIDTH = 5; // 5 bases on each side of the indel
	private int MISMATCH_CUTOFF = 1000000;
	private double AV_MISMATCHES_PER_READ = 1.5;
	
	
	private static String annGenomic = "GENOMIC";
	private static String annIntron = "INTRON";
	private static String annUTR = "UTR";
	private static String annCoding = "CODING";
	private static String annUnknown = "UNKNOWN";
	
	private SAMRecord lastRead;
	
	// "/humgen/gsa-scr1/GATK_Data/refGene.sorted.txt"
	
	@Override
	public void initialize() {
		coverage = new RunningCoverage(0,WINDOW_SIZE);
		
		if ( RefseqFileName != null ) {
			ReferenceOrderedData<rodRefSeq> refseq = new ReferenceOrderedData<rodRefSeq>("refseq",
					new java.io.File(RefseqFileName),rodRefSeq.class);
		
			refseqIterator = refseq.iterator();
			System.out.println("Using RefSeq annotations from "+RefseqFileName);
		}
		
		if ( refseqIterator == null ) System.out.println("No annotations available");
		
		int nSams = getToolkit().getArguments().samFiles.size(); 
		
		location = GenomeLocParser.createGenomeLoc(0,1);
		
		List<Set<String>> readGroupSets = getToolkit().getMergedReadGroupsByReaders();

		if ( call_somatic ) {
			if ( nSams != 2 ) {
				System.out.println("In --somatic mode two input bam files must be specified (normal/tumor)");
				System.exit(1);
			}
			normal_coverage = new RunningCoverage(0,WINDOW_SIZE);
			
			normalReadGroups = readGroupSets.get(0); // first -I option must specify normal.bam
			tumorReadGroups = readGroupSets.get(1); // second -I option must specify tumor.bam
		} else {
			if ( nSams != 1 ) System.out.println("WARNING: multiple input files specified. \n"+
					"WARNING: Without --somatic option they will be merged and processed as a single sample");
		}
/*
		List<Set<String>> sample_sets = getToolkit().getSamplesByReaders();
		for ( int i = 0 ; i < sample_sets.size() ; i++ ) {
			System.out.print("Reader "+i);
			for ( String s : sample_sets.get(i) ) System.out.print(" " + s);
			System.out.println();
		}

		List<Set<String>> lib_sets = getToolkit().getLibrariesByReaders();
		for ( int i = 0 ; i < lib_sets.size() ; i++ ) {
			System.out.print("Reader "+i);
			for ( String s : lib_sets.get(i) ) System.out.print(" " + s);
			System.out.println();
		}
*/
/*		
		for ( int i = 0 ; i < readGroupSets.size() ; i++ ) {
			System.out.print("Reader "+i);
			for ( String s : readGroupSets.get(i) ) System.out.print(" " + s);
			System.out.println();
		}
		assignReadGroups1();
*/
		try {
			output = new java.io.FileWriter(bed_file);
		} catch (IOException e) {
			throw new StingException("Failed to open file for writing BED output");
		}
	}

	/*
	void assignReadGroups(final SAMFileHeader mergedHeader) {
		
		Set<String> normal = new HashSet<String>(); // list normal samples here
		Set<String> tumor = new HashSet<String>(); // list tumor samples here
		
		SAMFileReader rn = new SAMFileReader(getToolkit().getArguments().samFiles.get(0));
		for ( SAMReadGroupRecord rec : rn.getFileHeader().getReadGroups() ) {
			normal.add(new String(rec.getSample()));
		}
		rn.close();
		rn = new SAMFileReader(getToolkit().getArguments().samFiles.get(1));
		for ( SAMReadGroupRecord rec : rn.getFileHeader().getReadGroups() ) {
			tumor.add(new String(rec.getSample()));
		}
		rn.close();		

		// now we know what samples are normal, and what are tumor; let's assign dynamic read groups we get in merged header:
		for ( SAMReadGroupRecord mr : mergedHeader.getReadGroups() ) {
			if ( normal.contains(mr.getSample() ) ) {
				normal_samples.add( new String(mr.getReadGroupId()) );
				System.out.println("Read group "+ mr.getReadGroupId() + "--> Sample "+ mr.getSample() + " (normal)");
			} else if ( tumor.contains(mr.getSample() ) ) {
				tumor_samples.add( new String(mr.getReadGroupId()) );
				System.out.println("Read group "+ mr.getReadGroupId() + "--> Sample "+ mr.getSample() + " (tumor)");
			} else throw new StingException("Unrecognized sample "+mr.getSample() +" in merged SAM stream");
		}
		System.out.println();
		
	}
	
*/
	
	
	@Override
	public Integer map(char[] ref, SAMRecord read) {
		

		if ( read.getReadUnmappedFlag() ||
			 read.getDuplicateReadFlag() ||
			 read.getNotPrimaryAlignmentFlag() ||
			 read.getMappingQuality() == 0 ) { 
			return 0; // we do not need those reads!
		}
		
		if ( read.getReferenceIndex() != currentContigIndex ) { 
			// we just jumped onto a new contig
			
			if ( read.getReferenceIndex() < currentContigIndex ) // paranoidal 
				throw new StingException("Read "+read.getReadName()+": contig is out of order");
			
			if ( call_somatic) emit_somatic(1000000000, true); // print remaining indels from the previous contig (if any);
			else emit(1000000000,true);
			currentContigIndex = read.getReferenceIndex();
			currentPosition = read.getAlignmentStart();
			refName = new String(read.getReferenceName());

			location = GenomeLocParser.setContig(location,refName);
			
			coverage.clear(); // reset coverage window; this will also set reference position to 0
			if ( call_somatic) normal_coverage.clear();
		}
		
		if ( read.getAlignmentStart() < currentPosition ) 
			throw new StingException("Read "+read.getReadName() +" out of order on the contig\n"+
					"Read starts at "+refName+":"+read.getAlignmentStart()+"; last read seen started at "+refName+":"+currentPosition
					+"\nLast read was: "+lastRead.getReadName()+" RG="+lastRead.getAttribute("RG")+" at "+lastRead.getAlignmentStart()+"-"
					+lastRead.getAlignmentEnd()+" cigar="+lastRead.getCigarString());
		
		currentPosition = read.getAlignmentStart();
		
		if ( read.getAlignmentStart() < coverage.getStart() ) {
			// should never happen
			throw new StingException("Read "+read.getReadName()+": out of order on the contig\n"+
					"Read starts at "+read.getReferenceName()+":"+read.getAlignmentStart()+ " (cigar="+read.getCigarString()+
					"); window starts at "+coverage.getStart());
		}
		
		lastRead = read;
		
		// a little trick here: we want to make sure that current read completely fits into the current
		// window so that we can accumulate the coverage/indel counts over the whole length of the read.
		// The ::getAlignmentEnd() method returns the last position on the reference where bases from the
		// read actually match (M or D cigar elements). After our cleaning procedure, we can have reads that end
		// with I element, which is not gonna be counted into alignment length on the reference. On the other hand, 
		// in this program we assign insertions, internally, to the first base *after* the insertion position. 
		// Hence, we have to make sure that that extra base is already in the window or we will get IndexOutOfBounds.
		
		long alignmentEnd = read.getAlignmentEnd();
		Cigar c = read.getCigar();
		if ( c.getCigarElement(c.numCigarElements()-1).getOperator() == CigarOperator.I) alignmentEnd++;
		
		if ( alignmentEnd > coverage.getStop()) {
			
			// we don't emit anything until we reach a read that does not fit into the current window.
			// At that point we shift the window to the start of that read and emit everything prior to 
			// that position (reads are sorted, so we are not gonna see any more coverage at those lower positions).
			// Clearly, we assume here that window is large enough to accomodate any single read, so simply shifting
			// the window to the read's start will ensure that the read fits...
			
			if ( call_somatic ) emit_somatic( read.getAlignmentStart(), false );
			else emit( read.getAlignmentStart(), false );

			if ( read.getAlignmentEnd() > coverage.getStop()) {
				// ooops, looks like the read does not fit into the current window!!
				throw new StingException("Read "+read.getReadName()+": out of coverage window bounds.Probably window is too small.\n"+
					"Read length="+read.getReadLength()+"; cigar="+read.getCigarString()+"; start="+
					read.getAlignmentStart()+"; end="+read.getAlignmentEnd()+"; window start="+coverage.getStart()+
					"; window end="+coverage.getStop());
			}
		}
		
		if ( call_somatic ) {
						
			String rg = (String)read.getAttribute("RG");
			if ( rg == null ) throw new StingException("Read "+read.getReadName()+" has no read group in merged stream");
			
			if ( normalReadGroups.contains(rg) ) {
				normal_coverage.add(read,ref);
			} else if ( tumorReadGroups.contains(rg) ) {
				coverage.add(read,ref);
			} else {
				throw new StingException("Unrecognized read group in merged stream: "+rg);
			}
		} else { 
			coverage.add(read, ref);
		}
		
		return 1;
	}

	/** Returns the indel variant with the largest count (ie consensus) among all the observed
	 * variants, and the total count of all observations of any indels (including non-consensus) 
	 * @param variants
	 * @return
	 */
	private Pair<IndelVariant,Integer> findConsensus(List<IndelVariant> variants) {
		int total_variant_count = 0;
		int max_variant_count = 0;
		IndelVariant v = null;
		
		for ( IndelVariant var : variants ) {
			int cnt = var.getCount();
			total_variant_count +=cnt;
			if ( cnt > max_variant_count ) v = var;
		}
		return new Pair<IndelVariant,Integer>(v,total_variant_count);
	}	
	
	/** Returns true if consensus (specified by the pair) should be considered a call given current values
	 * of the fraction cutoffs
	 * @param p
	 * @param coverage
	 * @return
	 */
	private boolean isCall(Pair<IndelVariant,Integer> p, int coverage) {
		return ( (double)p.first.getCount() > minFraction * coverage && (double) p.first.getCount() > minConsensusFraction*p.second ); 		
	}

	/** Build output line for bed file and write it to the specified output writer if the latter is not null; 
	 * the line is also returned by this method as a String
	 * 
	 * @param p
	 * @param coverage
	 * @param pos
	 * @param bedOutput
	 * @return
	 */
	private String makeBedLine(Pair<IndelVariant,Integer> p, int coverage, long pos, java.io.Writer bedOutput) {
		int event_length = p.first.lengthOnRef();
		if ( event_length < 0 ) event_length = 0;
		StringBuffer message = new StringBuffer();
        message.append(refName+"\t"+(pos-1)+"\t");
        if ( FORMAT_1KG )
            message.append(p.first.getBases().length() + "\t" + (event_length > 0 ? "D" : "I") + "\t" + p.first.getBases() + "\t" + p.first.getSamples());
        else
            message.append((pos-1+event_length)+"\t"+(event_length>0? "-":"+")+p.first.getBases() +":"+p.second+"/"+coverage);

		if ( bedOutput != null ) {
			try {
				bedOutput.write(message.toString()+"\n");
			} catch (IOException e) {
				System.out.println(e.getMessage());
				e.printStackTrace();
				throw new StingException("Error encountered while writing into output BED file");
			}
		}
		return message.toString();
	}

	/** Same as makeBedLine(Pair,int,long,Writer), but only builds and returns the line without writing it anywhere.
	 * 
	 * @param p
	 * @param coverage
	 * @param pos
	 * @return
	 */
	private String makeBedLine(Pair<IndelVariant,Integer> p, int coverage, long pos) {
		return makeBedLine(p, coverage, pos, null);
	}

	/** Output indel calls up to the specified position and shift the coverage array(s): after this method is executed
	 * first elements of the coverage arrays map onto 'position'
	 * 
	 * @param position
	 */
	private void emit(long position, boolean force) {
		
		long move_to = position; // we will shift to position move_to; it's initialized with 'position',
								 // but it may end up being smaller (delayed shift), if we have not
			                     // covered MISMATCH_WIDTH bases to the right of the last indel yet.
		
		// walk along the coverage window and emit indels up to the position we are trying ot shift the window to
		for ( long pos = coverage.getStart() ; pos < Math.min(position,coverage.getStop()+1) ; pos++ ) {
			
			List<IndelVariant> variants = coverage.indelsAt(pos);
			if ( variants.size() == 0 ) continue; // no indels at current position
			
			// if we are here, we got a variant

			int cov = coverage.coverageAt(pos);
			
			if ( cov < minCoverage ) continue; // low coverage
			
//			System.out.println("indel at "+pos);
			
			// region around the current indel we need to have covered in order to compute mismatch rate:
			long left = Math.max( pos-MISMATCH_WIDTH, coverage.getStart() );
			long right = pos+MISMATCH_WIDTH;
			
			// if right < position we are shifting to, we already have all the coverage within MISMATCH_WINDOW bases around the indel,
			// so we can proceed with counting mismatches and emitting the indel
			
			if ( right >= position && ! force) { 
				// we are not asked to force-shift, and there is more coverage around the current indel that we need to collect
					                            
				// we are not asked to force-shift, and there's still additional coverage to the right of current indel, so its too early to emit it;
				// instead we shift only up to current indel pos - MISMATCH_WIDTH, so that we could keep collecting that coverage
				move_to = left; 
//				System.out.println("right="+right+" requested="+position+" stopped at="+left);
				break;
			}

			if ( right > coverage.getStop() ) right = coverage.getStop(); // if indel is too close to the end of the window but we need to emit anyway 
			
			// count mismatches around the current indel, inside specified window (MISMATCH_WIDTH on each side):
			int total_mismatches = 0;
			for ( long k = left; k <= right ; k++ ) total_mismatches+=coverage.mismatchesAt(k);
			
			if ( total_mismatches > MISMATCH_CUTOFF || total_mismatches > ((double)cov)*AV_MISMATCHES_PER_READ) {
				out.println(refName+"\t"+(pos-1)+"\t"+
				"\tTOO DIRTY\t"+total_mismatches);
				continue; // too dirty
			}
			
			location = GenomeLocParser.setStart(location,pos); location = GenomeLocParser.setStop(location,pos); // retrieve annotation data
			rodRefSeq annotation = (refseqIterator == null ? null : refseqIterator.seekForward(location));

			Pair<IndelVariant,Integer> p = findConsensus(variants);
			if ( isCall(p,cov) ) { 	
				String message = makeBedLine(p,cov,pos,output);
				String annotationString = (refseqIterator == null ? "" : getAnnotationString(annotation));
				
				if ( verbose ) out.println(message + "\t"+ annotationString);
			} 
//			else { System.out.println("not a call: count="+p.first.count+" total_count="+p.second+" cov="+cov+
//					" minConsensusF="+((double)p.first.count)/p.second+" minF="+((double)p.first.count)/cov); }
//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
		}
//		System.out.println("Shifting to " + move_to+" ("+position+")");
		coverage.shift((int)(move_to - coverage.getStart() ) );
	}

	/** Output somatic indel calls up to the specified position and shift the coverage array(s): after this method is executed
	 * first elements of the coverage arrays map onto 'position'
	 * 
	 * @param position
	 */
	private void emit_somatic(long position, boolean force) {
		
		long move_to = position;
		
		for ( long pos = coverage.getStart() ; pos < Math.min(position,coverage.getStop()+1) ; pos++ ) {
			
			
			List<IndelVariant> tumor_variants = coverage.indelsAt(pos);
			List<IndelVariant> normal_variants = normal_coverage.indelsAt(pos);

			if ( tumor_variants.size() == 0 ) continue; // no indels in tumor

					
			int tumor_cov = coverage.coverageAt(pos);
			int normal_cov = normal_coverage.coverageAt(pos);
			
			if ( tumor_cov < minCoverage ) continue; // low coverage
			if ( normal_cov < minNormalCoverage ) continue; // low coverage
			
			long left = Math.max( pos-MISMATCH_WIDTH, coverage.getStart() );
			long right = pos+MISMATCH_WIDTH;
			
			if ( right >= position && ! force) { 
				// we are not asked to force-shift, and there is more coverage around the current indel that we need to collect
					                            
				// we are not asked to force-shift, and there's still additional coverage to the right of current indel, so its too early to emit it;
				// instead we shift only up to current indel pos - MISMATCH_WIDTH, so that we could keep collecting that coverage
				move_to = left; 
				break;
			}

			if ( right > coverage.getStop() ) right = coverage.getStop(); // if indel is too close to the end of the window but we need to emit anyway 

			// count mismatches around the current indel, inside specified window (MISMATCH_WIDTH on each side):
			int total_mismatches_normal = 0;
			int total_mismatches_tumor = 0;
			for ( long k = left; k <= right ; k++ ) {
				total_mismatches_tumor+=coverage.mismatchesAt(k);
				total_mismatches_normal+=normal_coverage.mismatchesAt(k);
			}
			
			if ( total_mismatches_normal > MISMATCH_CUTOFF || total_mismatches_normal > ((double)normal_cov)*AV_MISMATCHES_PER_READ) {
				out.println(refName+"\t"+(pos-1)+"\t"+
				"\tNORMAL TOO DIRTY\t"+total_mismatches_normal);
				continue; // too dirty
			}
			if ( total_mismatches_tumor > MISMATCH_CUTOFF || total_mismatches_tumor > ((double)tumor_cov)*AV_MISMATCHES_PER_READ) {
				out.println(refName+"\t"+(pos-1)+"\t"+
				"\tTUMOR TOO DIRTY\t"+total_mismatches_tumor);
				continue; // too dirty
			}
			location = GenomeLocParser.setStart(location,pos); location = GenomeLocParser.setStop(location,pos); // retrieve annotation data
			rodRefSeq annotation = (refseqIterator == null ? null : refseqIterator.seekForward(location));
			
			Pair<IndelVariant,Integer> p_tumor = findConsensus(tumor_variants);
			if ( isCall(p_tumor,tumor_cov) ) { 	
				String message = makeBedLine(p_tumor,tumor_cov,pos);
				String annotationString = (refseqIterator == null ? "" : getAnnotationString(annotation));
				
			
				if ( normal_variants.size() == 0 ) { 		

					try {
						output.write(message+"\n");
					} catch (IOException e) {
						System.out.println(e.getMessage());
						e.printStackTrace();
						throw new StingException("Error encountered while writing into output BED file");
					}
					message += "\tSOMATIC\t0/"+normal_cov;
				} else {
					Pair<IndelVariant,Integer> p_normal = findConsensus(tumor_variants);
					
					message += "\tGERMLINE\t"+p_normal.second+"/"+normal_cov;
				}
				if ( verbose ) out.println(message + "\t"+ annotationString);
			}
//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
		}
		
		coverage.shift((int)(move_to - coverage.getStart() ) );
		normal_coverage.shift((int)(move_to - normal_coverage.getStart() ) );
	}


	private String getAnnotationString(rodRefSeq ann) {
		if ( ann == null ) return annGenomic;
		else {
			StringBuilder b = new StringBuilder();
			if ( ann.isExon() ) {
				if ( ann.isCoding() ) b.append(annCoding);
				else b.append(annUTR);
			} else {
				if ( ann.isCoding() ) b.append(annIntron);
				else b.append(annUnknown);
			}
			b.append('\t');
			Iterator<Transcript> it = ann.getTranscripts().iterator();
			b.append(it.next().getGeneName()); // there is at least one transcript in the list, guaranteed
//			while ( it.hasNext() ) { // 
//				t.getGeneName()
//			}
			return b.toString();
		}
		
	}
	
	@Override
	public void onTraversalDone(Integer result) {
		if ( call_somatic ) emit_somatic(1000000000, true);
		else emit(1000000000,true); // emit everything we might have left
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
        private HashSet<String> samples = new HashSet<String>();

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
			else return 0;
		}
		
		public void increment() { count+=1; }
		
        public void addSample(String sample) {
            if ( sample != null )
                samples.add(sample);
        }

        public String getSamples() {
            StringBuffer sb = new StringBuffer();
            Iterator<String> i = samples.iterator();
            while ( i.hasNext() ) {
                sb.append(i.next());
                if ( i.hasNext() )
                    sb.append(",");
            }
            return sb.toString();
        }

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
		private CircularArray.Int mismatches;
		
		private static List<IndelVariant> emptyIndelList;
		
		static {
			emptyIndelList = new ArrayList<IndelVariant>();
		}
		
		public RunningCoverage(long start, int length) {
			this.start = start;
			coverageWindow = new CircularArray.Int(length);
			indels = new CircularArray< List<IndelVariant> >(length);
			mismatches = new CircularArray.Int(length);
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
		
		public int mismatchesAt(final long refPos) { return mismatches.get((int)(refPos-start)); }
		
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
			final String readBases = r.getReadString().toUpperCase();

			
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
	 //      int mm=0;
	        
	        for ( int i = 0 ; i < nCigarElems ; i++ ) {

	            final CigarElement ce = c.getCigarElement(i);
	            IndelVariant.Type type = null;
	            String bases = null;
	            int eventPosition = posOnRef;
	            
	            
	            switch(ce.getOperator()) {
	            case I:
	                    type = IndelVariant.Type.I; 
	                    bases = readBases.substring(posOnRead,posOnRead+ce.getLength()); 
	                    // will increment position on the read below, there's no 'break' statement yet...
	            case H:
	            case S:
	            		// here we also skip hard and soft-clipped bases on the read; according to SAM format specification, 
               			// alignment start position on the reference points to where the actually aligned 
    					// (not clipped) bases go, so we do not need to increment reference position here	                    
	            		posOnRead += ce.getLength();
	                    break;
	            case D: 
	            		type = IndelVariant.Type.D;
	            		bases = new String( ref, posOnRef, ce.getLength() );
	                    posOnRef += ce.getLength();
	                    break;
	            case M: for ( int k = 0; k < ce.getLength(); k++, posOnRef++, posOnRead++ ) {
	            		    if ( readBases.charAt(posOnRead) != Character.toUpperCase(ref[posOnRef]) ) { // mismatch!
	            		    	mismatches.increment(localStart+posOnRef, 1); //mm++;
	            		    }
	            		}
	            		break; // advance along the gapless block in the alignment
	            default :
	                throw new IllegalArgumentException("Unexpected operator in cigar string: "+ce.getOperator());
	            }

	            if ( type == null ) continue; // element was not an indel, go grab next element...
	            
	            // we got an indel if we are here...
	            if ( i == 0 ) logger.warn("Indel at the start of the read "+r.getReadName());
	            if ( i == nCigarElems - 1) logger.warn("Indel at the end of the read "+r.getReadName());

	            try {
	            	// note that here we will be assigning indels to the first deleted base or to the first
	            	// base after insertion, not to the last base before the event!
	            	updateCount(localStart+eventPosition, type, bases, r);
	            } catch (IndexOutOfBoundsException e) {
					System.out.println("Read "+r.getReadName()+": out of coverage window bounds.Probably window is too small.\n"+
							"Read length="+r.getReadLength()+"; cigar="+r.getCigarString()+"; start="+
							r.getAlignmentStart()+"; end="+r.getAlignmentEnd()+"; window start="+getStart()+
							"; window end="+getStop());
	            	throw e;
	            }
	        }
	            
//			System.out.println(r.getReadName()+"\t"+(r.getReadNegativeStrandFlag()?"RC":"FW")+"\t"+r.getCigarString()+"\t"+mm);
//			System.out.println(AlignmentUtils.alignmentToString(r.getCigar(), readBases, new String(ref), 0));
			
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
		private void updateCount(int pos, IndelVariant.Type type, String bases, SAMRecord r) {
            List<IndelVariant> indelsAtSite = indels.get(pos);
            if ( indelsAtSite == null ) {
            	indelsAtSite = new ArrayList<IndelVariant>();
            	indels.set(pos, indelsAtSite);
            }
            
            String sample = null;
            Object readGroupAttr = r.getAttribute("RG");
            if ( readGroupAttr != null ) {
                SAMReadGroupRecord readGroup = r.getHeader().getReadGroup(readGroupAttr.toString());
                if ( readGroup != null ) {
                    Object readSampleAttr = readGroup.getAttribute("SM");
                    if ( readSampleAttr != null )
                        sample = readSampleAttr.toString();
                }
            }

            boolean found = false;
            for ( IndelVariant v : indelsAtSite ) {
            	if ( ! v.equals(type, bases) ) continue;

                v.increment();
                v.addSample(sample);
            	found = true;
            	break;
            }
            
            if ( ! found ) {
                IndelVariant v = new IndelVariant(type, bases);
                v.addSample(sample);
                indelsAtSite.add(v);

            }
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
			mismatches.shiftData(offset);
		}
	}

}
