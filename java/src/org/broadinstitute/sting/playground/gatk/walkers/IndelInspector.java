package org.broadinstitute.sting.playground.gatk.walkers;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.indels.*;
import net.sf.samtools.*;

import edu.mit.broad.picard.reference.ReferenceSequence;

// Investigates indels called in the alignment data

@WalkerName("InspectIndels")
public class IndelInspector extends ReadWalker<Integer, Integer> {
    @Argument(fullName="IndelVerbosity", shortName="iverb", required=false, defaultValue="SILENT", doc="Verbosity level: SILENT, PILESUMMARY, ALIGNMENTS") public String VERBOSITY_LEVEL;
    @Argument(fullName="OutputNormal", shortName="outNorm", required=true, doc="Output file (sam or bam) for non-indel related reads and indel reads that were not improved") public String OUT1;
    @Argument(fullName="OutputIndel", shortName="outIndel", required=true, doc="Output file (sam or bam) for improved (realigned) indel related reads") public String OUT2;
    @Argument(fullName="ControlRun", shortName="paranoid", required=false, defaultValue="false", doc="If true, all reads that would be otherwise picked and processed by this tool will be saved, unmodified, into the OutputNormal file") public boolean CONTROL_RUN;
    @Argument(fullName="ErrorCheckMode", shortName="error", required=false, doc="Error counting mode: MM (mismatches only (from sam tags)), MC (mismatches only doing actual mismatch count on the fly (use this if tags are incorrectly set)), ERR (errors (arachne style: mm+gap lengths)), MG (count mismatches and gaps as one error each)") public String ERR_MODE;
    @Argument(fullName="MaxErrors", shortName="maxError", required=false, defaultValue="10000", doc="Maximum number of errors allowed (see ERR_MODE)") public Integer MAX_ERRS;
    @Argument(fullName="MaxReadLength", shortName="maxRead", required=false, defaultValue="-1", doc="Ignore reads that are longer than the specified cutoff (not a good way to do things but might be necessary because of performance issues)") public int MAX_READ_LENGTH;

    private int discarded_cigar_count = 0;
    private int discarded_long_read_count = 0;
    private int discarded_maxerr_count = 0;
    private int unmapped_read_count = 0;
    private int total_reads = 0;
    private IndelRecordPileCollector collector;
    private PileBuilder pileBuilder = null;
    private PassThroughWriter ptWriter;
    private String cur_contig = null;
    
    
 //   public boolean filter(LocusContext context, SAMRecord read) {
  //      // we only want aligned reads
  //      return !read.getReadUnmappedFlag();
 //   }

    public void initialize() {
     	
        if ( ERR_MODE != null && ! ERR_MODE.equals("MM") && ! ERR_MODE.equals("MG") && ! ERR_MODE.equals("ERR") && ! ERR_MODE.equals("MC") ) {
            err.println("Unknown value specified for ERR_MODE: "+ERR_MODE+ "... skipping it.");
            ERR_MODE = null;
        }

        SAMFileHeader header = getToolkit().getSamReader().getFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        ptWriter = new PassThroughWriter(OUT1, header);
        if ( ! CONTROL_RUN ) pileBuilder = new PileBuilder(OUT2, header, ptWriter);

        try {
            if ( CONTROL_RUN ) collector = new IndelRecordPileCollector(ptWriter, new DiscardingPileReceiver() );
            else collector = new IndelRecordPileCollector(ptWriter, pileBuilder );
        } catch(Exception e) {
            err.println(e.getMessage());
        }
        if ( collector == null )
            return;

        collector.setControlRun(CONTROL_RUN);

        if ( ! CONTROL_RUN ) {
            if ( VERBOSITY_LEVEL.toUpperCase().equals("SILENT")) pileBuilder.setVerbosity(PileBuilder.SILENT);
            else if ( VERBOSITY_LEVEL.toUpperCase().equals("PILESUMMARY") ) pileBuilder.setVerbosity(PileBuilder.PILESUMMARY);
            else if ( VERBOSITY_LEVEL.toUpperCase().equals("ALIGNMENTS") ) pileBuilder.setVerbosity(PileBuilder.ALIGNMENTS);
            else {
                err.println("Unrecognized VERBOSITY_LEVEL setting... skipping it.");
                pileBuilder.setVerbosity(PileBuilder.SILENT);
            }
        }
    }

    public Integer map(LocusContext context, SAMRecord read) {

    	total_reads++;
 
    	if ( read.getReadUnmappedFlag() ) {
    		unmapped_read_count++;
    		ptWriter.receive(read); // make sure we still write the read to the output, we do not want to lose data!
    		return 0;
    	}
    	
        ReferenceSequence contig_seq = context.getReferenceContig();
        if ( read.getReferenceName() != cur_contig) {
            cur_contig = read.getReferenceName();
            out.println("Contig "+cur_contig);
            String refstr = new String(contig_seq.getBases());
             if ( !CONTROL_RUN )
                pileBuilder.setReferenceSequence(refstr);
            out.println("loaded contig "+cur_contig+" (index="+read.getReferenceIndex()+"); length="+contig_seq.getBases().length+" tst="+contig_seq.toString());
        }

        //    if ( cur_contig.equals("chrM") || GenomeLoc.compareContigs(cur_contig,"chrY") > 0 ) continue; // skip chrM and unplaced contigs for now

        if ( MAX_READ_LENGTH > 0 && read.getReadLength() > MAX_READ_LENGTH ) {
        		ptWriter.receive(read);
                discarded_long_read_count++;
                return 0;
        }

        // we currently do not know how to deal with cigars containing elements other than M,I,D, so
        // let's just skip the reads that contain those other elements (clipped reads?)
        Cigar c = read.getCigar();
        boolean cigar_acceptable = true;
        for ( int z = 0 ; z < c.numCigarElements() ; z++ ) {
            CigarElement ce = c.getCigarElement(z);
            switch ( ce.getOperator() ) {
            case M:
            case I:
            case D: break;
            default:
                cigar_acceptable = false;
            }
        }
        if ( ! cigar_acceptable ) {
        		ptWriter.receive(read);
                discarded_cigar_count++;
                return 0;
        }

        int err = -1;
/*
            System.out.println("MM:     "+numMismatches(r));
            System.out.println("direct: "+numMismatchesDirect(r,contig_seq));
            System.out.print("  ");
            for ( int i = read.getAlignmentStart() - 1 ; i < read.getAlignmentEnd() ; i++ ) System.out.print((char)contig_seq.getBases()[i]);
            System.out.println();
            System.out.println((read.getReadNegativeStrandFlag()?"<-":"->")+read.getReadString());
            System.out.println("cigar: "+read.getCigarString());
            System.out.println();
            if (counter++ == 20 ) break;
            continue;
*/

        if ( ERR_MODE != null ) {
            if ( ERR_MODE.equals("MM"))       err = numMismatches(read,contig_seq);
            else if ( ERR_MODE.equals("MC") ) err = AlignmentUtils.numMismatches(read,contig_seq);
            else if ( ERR_MODE.equals("ERR")) err = numErrors(read,contig_seq);
            else if ( ERR_MODE.equals("MG"))  err = numMismatchesGaps(read,contig_seq);
            if ( err > MAX_ERRS ) {
            	ptWriter.receive(read);
            	discarded_maxerr_count++;
                return 0;
            }
        }

        collector.receive(read);
        return 1;
    }

	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD
	 *
	 * @param r SAM record that must specify an alignment
	 * @return number of errors (number of mismatches plus total length of all insertions/deletions
	 * @throws RuntimeException
	 */
    private static int numErrors(SAMRecord r, ReferenceSequence refseq) throws RuntimeException {

		// NM currently stores the total number of mismatches in all blocks + 1
        int errs = numMismatches(r,refseq);

		// now we have to add the total length of all indels:
		Cigar c = r.getCigar();
		for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
			CigarElement ce = c.getCigarElement(i);
			switch( ce.getOperator()) {
			case M : break; // we already have correct number of mismatches
			case I :
			case D :
					errs += ce.getLength();
					break;
			default: throw new RuntimeException("Unrecognized cigar element");
			}
		}
		return errs;
	}

	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD
	 *
	 * @param r SAM record that must specify an alignment
	 * @return number of errors (number of mismatches plus total number of all insertions/deletions (each insertion or
	 * deletion will be counted as a single error regardless of the length)
	 * @throws RuntimeException
	 */
    private static int numMismatchesGaps(SAMRecord r,ReferenceSequence refseq) throws RuntimeException {

		// NM currently stores the total number of mismatches in all blocks + 1
        int errs = numMismatches(r,refseq);

		// now we have to add the total length of all indels:
		Cigar c = r.getCigar();
		for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
			CigarElement ce = c.getCigarElement(i);
			switch( ce.getOperator()) {
			case M : break; // we already have correct number of mismatches
			case I :
			case D :
					errs++;
					break;
			default: throw new RuntimeException("Unrecognized cigar element");
			}
		}
		return errs;
	}

    /** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD */
    private static int numMismatches(SAMRecord r, ReferenceSequence refseq) throws RuntimeException {

        // NM currently stores the total number of mismatches in all blocks + 1
        Integer i = (Integer)r.getAttribute("NM");
        if ( i == null ) return AlignmentUtils.numMismatches(r,refseq);
        return ((Integer)r.getAttribute("NM")).intValue() - 1;

    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer result) {
    	collector.close(); // write the reads collector might be still holding
    	
        if ( ! CONTROL_RUN ) {
            pileBuilder.printStats();
            pileBuilder.close();
        }
        out.println("done.");
        out.println("Total reads processed: "+ total_reads);
        out.println("Unmodified: "+ ptWriter.getNumReadsReceived());
        if ( pileBuilder != null) out.println("Written into realigned piles: "+ pileBuilder.getNumReadsWritten());
        if ( pileBuilder != null) out.println("Received for realignment: "+ pileBuilder.getNumReadsReceived());
        out.println("Discarded reads with non-M,I,D cigar elements: "+ discarded_cigar_count);
        out.println("Discarded long reads (above "+MAX_READ_LENGTH+" bp): "+ discarded_long_read_count);
        out.println("Discarded reads with error counts exceeding the threshold: "+ discarded_maxerr_count);
        out.println("Unmapped reads (passed through to the output): "+ unmapped_read_count);
        out.println();
        collector.printLengthHistograms();
        ptWriter.close();    
    }
}