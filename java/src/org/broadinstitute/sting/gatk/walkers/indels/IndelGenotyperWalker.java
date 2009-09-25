package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.*;

import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.filters.PlatformUnitFilter;
import org.broadinstitute.sting.gatk.filters.PlatformUnitFilterHelper;
import org.broadinstitute.sting.playground.utils.CircularArray;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;


@ReadFilters({Platform454Filter.class, ZeroMappingQualityReadFilter.class, PlatformUnitFilter.class})
public class IndelGenotyperWalker extends ReadWalker<Integer,Integer> {
    @Argument(fullName="outputFile", shortName="O", doc="output file name (defaults to BED format)", required=true)
    java.io.File bed_file;
    @Argument(fullName="1kg_format", shortName="1kg", doc="output in 1000 genomes format", required=false)
    boolean FORMAT_1KG = false;
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
			doc="Minimum fraction of reads with CONSENSUS indel at a site, out of all reads covering the site, required for a consensus call"+
			" (fraction of non-consensus indels at the site is not considered here, see minConsensusFraction)", required=false)
	double minFraction = 0.3;
	@Argument(fullName="minConsensusFraction", shortName="minConsensusFraction", 
			doc="Minimum fraction of CONSENSUS indel observations at a site wrt all indel observations at the site required to make the call", required=false)
	double minConsensusFraction = 0.7;
	@Argument(fullName="minIndelCount", shortName="minCnt", 
			doc="Minimum count of reads supporting consensus indel required for making the call. "+
			" This filter supercedes minFraction, i.e. indels with acceptable minFraction at low coverage "+
			"(minIndelCount not met) will not pass.", required=false)
	int minIndelCount = 0;
	@Argument(fullName="refseq", shortName="refseq", 
			doc="Name of RefSeq transcript annotation file. If specified, indels will be annotated as GENOMIC/UTR/INTRON/CODING", required=false)
	String RefseqFileName = null;
    @Argument(fullName="blacklistedLanes", shortName="BL",
            doc="Name of lanes (platform units) that should be ignored. Reads coming from these lanes will never be seen "+
                    "by this application, so they will not contribute indels to consider and will not be counted.", required=false)
        PlatformUnitFilterHelper dummy;
     @Argument(fullName="indel_debug", shortName="idebug", doc="Detailed printout for debugging",required=false) Boolean DEBUG = false;
	 
	private static int WINDOW_SIZE = 200;
	private RunningCoverage tumor_coverage;
	private RunningCoverage normal_coverage; // when performing somatic calls, we will be using this one for normal, and 'tumor_coverage' for tumor
	private int currentContigIndex = -1;
	private int currentPosition = -1; // position of the last read we've seen on the current contig
	private String refName = null;
	private java.io.Writer output = null;
	private GenomeLoc location = null;
	
    private SeekableRODIterator<rodRefSeq> refseqIterator=null;

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
		normal_coverage = new RunningCoverage(0,WINDOW_SIZE);
		
		if ( RefseqFileName != null ) {
			ReferenceOrderedData<rodRefSeq> refseq = new ReferenceOrderedData<rodRefSeq>("refseq",
					new java.io.File(RefseqFileName), rodRefSeq.class);
		
			refseqIterator = refseq.iterator();
			logger.info("Using RefSeq annotations from "+RefseqFileName);
		}
		
		if ( refseqIterator == null ) logger.info("No annotations available");
		
		int nSams = getToolkit().getArguments().samFiles.size(); 
		
		location = GenomeLocParser.createGenomeLoc(0,1);
		
		List<Set<String>> readGroupSets = getToolkit().getMergedReadGroupsByReaders();

		if ( call_somatic ) {
			if ( nSams != 2 ) {
				System.out.println("In --somatic mode two input bam files must be specified (normal/tumor)");
				System.exit(1);
			}
			tumor_coverage = new RunningCoverage(0,WINDOW_SIZE);
			
			normalReadGroups = readGroupSets.get(0); // first -I option must specify normal.bam
			tumorReadGroups = readGroupSets.get(1); // second -I option must specify tumor.bam
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

        if ( DEBUG ) {
//            System.out.println("DEBUG>> read at "+ read.getAlignmentStart()+"-"+read.getAlignmentEnd()+
//                    "("+read.getCigarString()+")");
            if ( read.getDuplicateReadFlag() ) System.out.println("DEBUG>> Duplicated read (IGNORED)");
        }

		if ( AlignmentUtils.isReadUnmapped(read) ||
			 read.getDuplicateReadFlag() ||
			 read.getNotPrimaryAlignmentFlag() ||
			 read.getMappingQuality() == 0 ) { 
			return 0; // we do not need those reads!
		}
		
		if ( read.getReferenceIndex() != currentContigIndex ) { 
			// we just jumped onto a new contig
			
			if ( read.getReferenceIndex() < currentContigIndex ) // paranoidal 
				throw new StingException("Read "+read.getReadName()+": contig is out of order; input BAM file is unsorted");
			
            // print remaining indels from the previous contig (if any);
 			if ( call_somatic ) emit_somatic(1000000000, true);
			else emit(1000000000,true);

			currentContigIndex = read.getReferenceIndex();
			currentPosition = read.getAlignmentStart();
			refName = new String(read.getReferenceName());

			location = GenomeLocParser.setContig(location,refName);
			
			normal_coverage.clear(); // reset coverage window; this will also set reference position to 0
			if ( call_somatic) tumor_coverage.clear();
		}

        // we have reset the window to the new contig if it was required and emitted everything we collected
        // on a previous contig. At this point we are guaranteed that we are set up properly for working
        // with the contig of the current read.

        // NOTE: all the sanity checks and error messages below use normal_coverage only. We make sure that normal_coverage and
        // tumor_coverage are synchronized exactly (windows are always shifted together by emit_somatic), so it's safe

		if ( read.getAlignmentStart() < currentPosition ) // oops, read out of order?
			throw new StingException("Read "+read.getReadName() +" out of order on the contig\n"+
					"Read starts at "+refName+":"+read.getAlignmentStart()+"; last read seen started at "+refName+":"+currentPosition
					+"\nLast read was: "+lastRead.getReadName()+" RG="+lastRead.getAttribute("RG")+" at "+lastRead.getAlignmentStart()+"-"
					+lastRead.getAlignmentEnd()+" cigar="+lastRead.getCigarString());
		
		currentPosition = read.getAlignmentStart();
		
		if ( read.getAlignmentStart() < normal_coverage.getStart() ) {
			// should never happen
			throw new StingException("Read "+read.getReadName()+": out of order on the contig\n"+
					"Read starts at "+read.getReferenceName()+":"+read.getAlignmentStart()+ " (cigar="+read.getCigarString()+
					"); window starts at "+normal_coverage.getStart());
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
		
		if ( alignmentEnd > normal_coverage.getStop()) {
			
			// we don't emit anything until we reach a read that does not fit into the current window.
			// At that point we shift the window to the start of that read and emit everything prior to 
			// that position (reads are sorted, so we are not gonna see any more coverage at those lower positions).
			// Clearly, we assume here that window is large enough to accomodate any single read, so simply shifting
			// the window to the read's start will ensure that the read fits...
			
			if ( call_somatic ) emit_somatic( read.getAlignmentStart(), false );
			else emit( read.getAlignmentStart(), false );

			if ( read.getAlignmentEnd() > normal_coverage.getStop()) {
				// ooops, looks like the read does not fit into the window even after the latter was shifted!!
				throw new StingException("Read "+read.getReadName()+": out of coverage window bounds. Probably window is too small.\n"+
					"Read length="+read.getReadLength()+"; cigar="+read.getCigarString()+"; start="+
					read.getAlignmentStart()+"; end="+read.getAlignmentEnd()+"; window start (after trying to accomodate the read)="+normal_coverage.getStart()+
					"; window end="+normal_coverage.getStop());
			}
		}
		
		if ( call_somatic ) {
						
			String rg = (String)read.getAttribute("RG");
			if ( rg == null ) throw new StingException("Read "+read.getReadName()+" has no read group in merged stream. RG is required for somatic calls.");
			
			if ( normalReadGroups.contains(rg) ) {
				normal_coverage.add(read,ref);
			} else if ( tumorReadGroups.contains(rg) ) {
				tumor_coverage.add(read,ref);
			} else {
				throw new StingException("Unrecognized read group in merged stream: "+rg);
			}
		} else { 
			normal_coverage.add(read, ref);
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
            if ( DEBUG ) System.out.println("DEBUG>> Variant "+var.getBases()+" (cnt="+var.getCount()+")");
			int cnt = var.getCount();
			total_variant_count +=cnt;
			if ( cnt > max_variant_count ) {
                v = var;
                max_variant_count = cnt;
            }
		}
        if ( DEBUG ) System.out.println("DEBUG>> Returning: "+v.getBases()+" (cnt="+v.getCount()+") with total count of "+total_variant_count);
		return new Pair<IndelVariant,Integer>(v,total_variant_count);
	}	
	
	/** Returns true if consensus (specified by the pair) should be considered a call given current values
	 * of the cutoffs.
	 * @param p pair with first element being the consensus indel variant, the second element being the <i>total</i> (consensus+others)
	 * 		count of indels at the site.
	 * @param coverage total coverage (number of spanning reads, including those with indel(s)) at the site.
	 * @return
	 */
	private boolean isCall(Pair<IndelVariant,Integer> p, int coverage) {
        boolean ret = ( p.first.getCount() >= minIndelCount &&
				(double)p.first.getCount() > minFraction * coverage &&
				(double) p.first.getCount() > minConsensusFraction*p.second );
        if ( DEBUG && ! ret ) System.out.println("DEBUG>>  NOT a call: count="+p.first.count+" total_count="+p.second+" cov="+coverage+
            " minConsensusF="+((double)p.first.count)/p.second+" minF="+((double)p.first.count)/coverage);
		return ret;
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

	/** Output indel calls up to the specified position and shift the coverage array: after this method is executed, the
	 * first element of the coverage array maps onto 'position', or a few bases to the left of 'position' if we may need more
     * reads to get full NQS-style statistics for an indel in the close proximity of 'position'.
	 * 
	 * @param position
	 */
	private void emit(long position, boolean force) {
		
		long move_to = position; // we will shift to position move_to; it's initialized with 'position',
								 // but it may end up being smaller (delayed shift), if we have not
			                     // covered MISMATCH_WIDTH bases to the right of the last indel yet.
		
//		boolean debug = false;
//		if ( coverage.getStart() <= 19661504 && coverage.getStop() >= 19661504 ) debug = true;
		
		if ( DEBUG ) System.out.println("DEBUG>> Window: ["+normal_coverage.getStart()+", "+normal_coverage.getStop()+"]; shift requested: to "+position);
		
		// walk along the coverage window and emit indels up to the position we are trying ot shift the window to
		for ( long pos = normal_coverage.getStart() ; pos < Math.min(position,normal_coverage.getStop()+1) ; pos++ ) {
			
			List<IndelVariant> variants = normal_coverage.indelsAt(pos);
			if ( variants.size() == 0 ) continue; // no indels at current position, go check next one
			
			// if we are here, we got a variant

			int cov = normal_coverage.coverageAt(pos); // depth of coverage
			
			if ( cov < minCoverage ) continue; // coverage too low to make a call
			
			// region around the current indel we need to have covered in order to compute mismatch rate:
			long left = Math.max( pos-MISMATCH_WIDTH, normal_coverage.getStart() );
			long right = pos+MISMATCH_WIDTH;
			
			if ( DEBUG ) System.out.println("DEBUG>> Indel at "+pos);
			
			if ( right >= position && ! force) {
				// we are not asked to force-shift, and there's still additional coverage to the right of current indel, so its too early to emit it;
				// instead we shift only up to current indel pos - MISMATCH_WIDTH, so that we could keep collecting that coverage
				move_to = left; 
				if ( DEBUG ) System.out.println("DEBUG>> waiting for coverage; actual shift performed to "+ left);
				break; // abort, don't output current indel - yet
			}

            // ok, right < position we are shifting to (or we force-shift), so we already have all the coverage within
            // MISMATCH_WINDOW bases around the indel;
            // we can proceed with counting mismatches and emitting the indel:

			if ( right > normal_coverage.getStop() ) right = normal_coverage.getStop(); // in case indel is too close to the end of the window but we need to emit (force-shift)
			
			// count mismatches around the current indel, inside the specified window (MISMATCH_WIDTH on each side):
			int total_mismatches = 0;
			for ( long k = left; k <= right ; k++ ) total_mismatches+=normal_coverage.mismatchesAt(k);
			
			if ( total_mismatches > MISMATCH_CUTOFF || total_mismatches > ((double)cov)*AV_MISMATCHES_PER_READ) {
				out.println(refName+"\t"+(pos-1)+"\t"+
				"\tTOO DIRTY\t"+total_mismatches);
				normal_coverage.indelsAt(pos).clear(); // we dealt with this indel; don't want to see it again
					                            // (we might otherwise in the case when 1) there is another indel that follows
												// within MISMATCH_WIDTH bases and 2) we'd need to wait for more coverage for that next indel)
					                            
				continue; // too dirty
			}
			
			location = GenomeLocParser.setStart(location,pos); location = GenomeLocParser.setStop(location,pos); // retrieve annotation data
			RODRecordList<rodRefSeq> annotationList = (refseqIterator == null ? null : refseqIterator.seekForward(location));

			Pair<IndelVariant,Integer> p = findConsensus(variants);
			if ( isCall(p,cov) ) { 	
				String message = makeBedLine(p,cov,pos,output);
				String annotationString = (refseqIterator == null ? "" : getAnnotationString(annotationList));
				
				if ( verbose ) out.println(message + "\t"+ annotationString);
			}
			normal_coverage.indelsAt(pos).clear(); // we dealt with this indel; don't want to see it again 
                                            // (we might otherwise in the case when 1) there will be another indel that follows
			                                // within MISMATCH_WIDTH bases and 2) we'd need to wait for more coverage for that next indel)

//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
		}
		if ( DEBUG ) System.out.println("DEBUG>> Actual shift to " + move_to+" ("+position+")");
		normal_coverage.shift((int)(move_to - normal_coverage.getStart() ) );
	}

	/** Output somatic indel calls up to the specified position and shift the coverage array(s): after this method is executed
	 * first elements of the coverage arrays map onto 'position', or a few bases prior to the specified position
     * if there is an indel in close proximity to 'position' so that we may get more coverage around it later.
	 * 
	 * @param position
	 */
	private void emit_somatic(long position, boolean force) {
		
		long move_to = position;
		
		for ( long pos = tumor_coverage.getStart() ; pos < Math.min(position,tumor_coverage.getStop()+1) ; pos++ ) {
			
			
			List<IndelVariant> tumor_variants = tumor_coverage.indelsAt(pos);
			List<IndelVariant> normal_variants = normal_coverage.indelsAt(pos);

			if ( tumor_variants.size() == 0 ) continue; // no indels in tumor

					
			int tumor_cov = tumor_coverage.coverageAt(pos);
			int normal_cov = normal_coverage.coverageAt(pos);
			
			if ( tumor_cov < minCoverage ) {
                if ( DEBUG ) {
                    System.out.println("DEBUG>> Indel in tumor at "+pos+"; coverare in tumor="+tumor_cov+" (SKIPPED)");
                }
                continue; // low coverage
            }
			if ( normal_cov < minNormalCoverage ) {
                if ( DEBUG ) {
                    System.out.println("DEBUG>> Indel in tumor at "+pos+"; coverare in normal="+normal_cov+" (SKIPPED)");
                }
                continue; // low coverage
            }
			
            if ( DEBUG ) System.out.println("DEBUG>> Indel in tumor at "+pos);

			long left = Math.max( pos-MISMATCH_WIDTH, tumor_coverage.getStart() );
			long right = pos+MISMATCH_WIDTH;
			
			if ( right >= position && ! force) { 
				// we are not asked to force-shift, and there is more coverage around the current indel that we still need to collect
					                            
				// we are not asked to force-shift, and there's still additional coverage to the right of current indel, so its too early to emit it;
				// instead we shift only up to current indel pos - MISMATCH_WIDTH, so that we could keep collecting that coverage
				move_to = left; 
                if ( DEBUG ) System.out.println("DEBUG>> waiting for coverage; actual shift performed to "+ left);
				break;
			}

			if ( right > tumor_coverage.getStop() ) right = tumor_coverage.getStop(); // if indel is too close to the end of the window but we need to emit anyway (force-shift), adjust right

			// count mismatches around the current indel, inside specified window (MISMATCH_WIDTH on each side):
			int total_mismatches_normal = 0;
			int total_mismatches_tumor = 0;
			for ( long k = left; k <= right ; k++ ) {
				total_mismatches_tumor+=tumor_coverage.mismatchesAt(k);
				total_mismatches_normal+=normal_coverage.mismatchesAt(k);
			}
			
			if ( total_mismatches_normal > MISMATCH_CUTOFF || total_mismatches_normal > ((double)normal_cov)*AV_MISMATCHES_PER_READ) {
				out.println(refName+"\t"+(pos-1)+"\t"+
				"\tNORMAL TOO DIRTY\t"+total_mismatches_normal);
				tumor_coverage.indelsAt(pos).clear();
				normal_coverage.indelsAt(pos).clear(); 
					// we dealt with this indel; don't want to see it again 
                	// (we might otherwise in the case when 1) there is another indel that follows
					// within MISMATCH_WIDTH bases and 2) we'd need to wait for more coverage for that next indel)
				continue; // too dirty
			}
			if ( total_mismatches_tumor > MISMATCH_CUTOFF || total_mismatches_tumor > ((double)tumor_cov)*AV_MISMATCHES_PER_READ) {
				out.println(refName+"\t"+(pos-1)+"\t"+
				"\tTUMOR TOO DIRTY\t"+total_mismatches_tumor);
				tumor_coverage.indelsAt(pos).clear();
				normal_coverage.indelsAt(pos).clear(); 
					// we dealt with this indel; don't want to see it again 
                	// (we might otherwise in the case when 1) there is another indel that follows
					// within MISMATCH_WIDTH bases and 2) we'd need to wait for more coverage for that next indel)
				continue; // too dirty
			}
			location = GenomeLocParser.setStart(location,pos); location = GenomeLocParser.setStop(location,pos); // retrieve annotation data
            RODRecordList<rodRefSeq> annotationList = (refseqIterator == null ? null : refseqIterator.seekForward(location));
			
			Pair<IndelVariant,Integer> p_tumor = findConsensus(tumor_variants);
			if ( isCall(p_tumor,tumor_cov) ) { 	
				String message = makeBedLine(p_tumor,tumor_cov,pos);
				String annotationString = (refseqIterator == null ? "" : getAnnotationString(annotationList));
				
			
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
					Pair<IndelVariant,Integer> p_normal = findConsensus(normal_variants);
					
					message += "\tGERMLINE\t"+p_normal.second+"/"+normal_cov;
				}
				if ( verbose ) {
                    if ( refseqIterator == null ) out.println(message + "\t");
                    else out.println(message + "\t"+ annotationString);
                }
			}
			
			tumor_coverage.indelsAt(pos).clear();
			normal_coverage.indelsAt(pos).clear(); 
				// we dealt with this indel; don't want to see it again 
            	// (we might otherwise in the case when 1) there is another indel that follows
				// within MISMATCH_WIDTH bases and 2) we'd need to wait for more coverage for that next indel)
			
//			for ( IndelVariant var : variants ) {
//				System.out.print("\t"+var.getType()+"\t"+var.getBases()+"\t"+var.getCount());
//			}
		}
		
        if ( DEBUG ) System.out.println("DEBUG>> Actual shift to " + move_to+" ("+position+")");
		tumor_coverage.shift((int)(move_to - tumor_coverage.getStart() ) );
		normal_coverage.shift((int)(move_to - normal_coverage.getStart() ) );
	}


	private String getAnnotationString(RODRecordList<rodRefSeq> ann) {
		if ( ann == null ) return annGenomic;
		else {
			StringBuilder b = new StringBuilder();

			if ( rodRefSeq.isExon(ann) ) {
				if ( rodRefSeq.isCoding(ann) ) b.append(annCoding); // both exon and coding = coding exon sequence
				else b.append(annUTR); // exon but not coding = UTR
			} else {
				if ( rodRefSeq.isCoding(ann) ) b.append(annIntron); // not in exon, but within the coding region = intron
				else b.append(annUnknown); // we have no idea what this is. this may actually happen when we have a fully non-coding exon...
			}
			b.append('\t');
			b.append(((Transcript)ann.getRecords().get(0)).getGeneName()); // there is at least one transcript in the list, guaranteed
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
        public void increment() { count+=1; }
		
		/** Returns length of the event on the reference (number of deleted bases
		 * for deletions, -1 for insertions. 
		 * @return
		 */
		public int lengthOnRef() { 
			if ( type == Type.D ) return bases.length();
			else return 0;
		}
		

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

        // Lists will exactly mimic the reads covering corresponding base, in the right order;
        // value = 1 if read has a mismatch, 0 otherwise
        private CircularArray< List<Integer> > mm_flags;

        // Lists will exactly mimic the reads covering corresponding base, in the right order;
        // i-th value = base quality at this location in the i-th read 
        private CircularArray< List<Integer> > base_quals;
		
		private static List<IndelVariant> emptyIndelList;
        private static Integer ZERO = new Integer(0);
        private static Integer ONE = new Integer(1);

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
