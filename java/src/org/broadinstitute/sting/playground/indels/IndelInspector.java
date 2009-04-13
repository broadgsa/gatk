package org.broadinstitute.sting.playground.indels;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.HashMap;


import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.reference.ReferenceSequenceFileWalker;
import edu.mit.broad.picard.reference.ReferenceSequence;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.GenomeLoc;

public class IndelInspector extends CommandLineProgram {

    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Investigates indels called in the alignment data\n";
    @Option(shortName="I", doc="SAM or BAM file for calling",optional=true) public File INPUT_FILE;
    @Option(shortName="L",doc="Genomic interval to run on, as contig[:start[-stop]]; whole genome if not specified", optional=true) public String GENOME_LOCATION;
    @Option(shortName="V",doc="Verbosity level: SILENT, PILESUMMARY, ALIGNMENTS", optional=true) public String VERBOSITY_LEVEL;
    @Option(doc="Output file (sam or bam) for non-indel related reads and indel reads that were not improved") public String OUT1; 
    @Option(doc="Output file (sam or bam) for improved (realigned) indel related reads") public String OUT2;
    @Option(doc="[paranoid] If true, all reads that would be otherwise picked and processed by this tool will be saved, unmodified, into OUT1", optional=true) public Boolean CONTROL_RUN;
    @Option(doc="Error counting mode: MM - mismatches only (from sam tags), MC - mismatches only doing actual mismatch count on the fly (use this if tags are incorrectly set); ERR - errors (arachne style: mm+gap lengths), MG - count mismatches and gaps as one error each") public String ERR_MODE;
    @Option(doc="Maximum number of errors allowed (see ERR_MODE)") public Integer MAX_ERRS;
    @Option(shortName="R", doc="Reference fasta or fasta.gz file") public File REF_FILE;
    @Option(doc="Ignore reads that are longer than the specified cutoff (not a good way to do things but might be necessary because of performance issues)", optional=true) public Integer MAX_READ_LENGTH;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new IndelInspector().instanceMain(argv));
    }
    
    protected int doWork() {

        int discarded_cigar_count = 0;
        int discarded_long_read_count = 0;

        ReferenceSequenceFileWalker reference = new ReferenceSequenceFileWalker(
                    REF_FILE
            );

        if ( reference.getSequenceDictionary() == null ) {
            System.out.println("No reference sequence dictionary found. Abort.");
        }

        GenomeLoc.setupRefContigOrdering(reference.getSequenceDictionary());
        GenomeLoc location = null;
        if ( GENOME_LOCATION != null ) {
            location = GenomeLoc.parseGenomeLoc(GENOME_LOCATION);
        }
        
        if ( ! ERR_MODE.equals("MM") && ! ERR_MODE.equals("MG") && ! ERR_MODE.equals("ERR") && ! ERR_MODE.equals("MC") ) {
            System.out.println("Unknown value specified for ERR_MODE: "+ERR_MODE);
            return 1;
        }

        final SAMFileReader samReader = new SAMFileReader(getInputFile(INPUT_FILE,"/broad/1KG/"));
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        //        setContigOrdering(samReader);


        if ( MAX_READ_LENGTH == null ) MAX_READ_LENGTH = 1000000000;
         
        ReferenceSequence contig_seq = null;

        IndelRecordPileCollector col = null;
        PassThroughWriter ptWriter = new PassThroughWriter(OUT1,samReader.getFileHeader());
        PileBuilder pileBuilder = null;
        if ( CONTROL_RUN == null ) CONTROL_RUN=false;
        if ( ! CONTROL_RUN ) pileBuilder = new PileBuilder(OUT2,samReader.getFileHeader(),ptWriter);

        try {
            if ( CONTROL_RUN ) col = new IndelRecordPileCollector(ptWriter, new DiscardingPileReceiver() );
            else col = new IndelRecordPileCollector(ptWriter, pileBuilder );
        } catch(Exception e) { System.err.println(e.getMessage()); }
        if ( col == null ) return 1; 

        col.setControlRun(CONTROL_RUN);

        if ( ! CONTROL_RUN ) {
            if ( VERBOSITY_LEVEL == null ) VERBOSITY_LEVEL = new String("SILENT");
            if ( VERBOSITY_LEVEL.toUpperCase().equals("SILENT")) pileBuilder.setVerbosity(PileBuilder.SILENT);
            else if ( VERBOSITY_LEVEL.toUpperCase().equals("PILESUMMARY") ) pileBuilder.setVerbosity(PileBuilder.PILESUMMARY);
            else if ( VERBOSITY_LEVEL.toUpperCase().equals("ALIGNMENTS") ) pileBuilder.setVerbosity(PileBuilder.ALIGNMENTS);
            else {
                System.out.println("Unrecognized VERBOSITY_LEVEL setting.");
                return 1;
            }
        }

        String cur_contig = null;
        int counter = 0;

        for ( SAMRecord r : samReader ) {

            if ( r.getReadUnmappedFlag() ) continue; 
            if ( r.getReferenceName() != cur_contig) {
                cur_contig = r.getReferenceName();
                System.out.println("Contig "+cur_contig);
                // if contig is specified and we are past that contig, we are done:
                if ( location != null && GenomeLoc.compareContigs(cur_contig, location.getContig()) == 1 ) break;
                if ( location == null || GenomeLoc.compareContigs(cur_contig, location.getContig()) == 0 ) {
                    contig_seq = reference.get(r.getReferenceIndex());
                    String refstr = new String(contig_seq.getBases());
                    col.setReferenceSequence(refstr);
                    if (!CONTROL_RUN) pileBuilder.setReferenceSequence(refstr);
                    System.out.println("loaded contig "+cur_contig+" (index="+r.getReferenceIndex()+"); length="+contig_seq.getBases().length+" tst="+contig_seq.toString());
                }
            }

            // if contig is specified and we did not reach it yet, skip the records until we reach that contig:
            if ( location != null && GenomeLoc.compareContigs(cur_contig, location.getContig()) == -1 ) continue;


            if ( location != null && r.getAlignmentEnd() < location.getStart() ) continue;

            // if stop position is specified and we are past that, stop reading:
            if ( location != null && r.getAlignmentStart() > location.getStop() ) break;

            //    if ( cur_contig.equals("chrM") || GenomeLoc.compareContigs(cur_contig,"chrY") > 0 ) continue; // skip chrM and unplaced contigs for now

            // we currently do not know how to deal with cigars containing elements other than M,I,D, so 
            // let's just skip the reads that contain those other elements (clipped reads?)
            Cigar c = r.getCigar();
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
               	discarded_cigar_count++;
               	continue;
            }
            
            if ( r.getReadLength() > MAX_READ_LENGTH ) {
            		discarded_long_read_count++;
            		continue;
            }

            int err = -1;
/*
            System.out.println("MM:     "+numMismatches(r));
            System.out.println("direct: "+numMismatchesDirect(r,contig_seq));
            System.out.print("  ");
            for ( int i = r.getAlignmentStart() - 1 ; i < r.getAlignmentEnd() ; i++ ) System.out.print((char)contig_seq.getBases()[i]);
            System.out.println();
            System.out.println((r.getReadNegativeStrandFlag()?"<-":"->")+r.getReadString());
            System.out.println("cigar: "+r.getCigarString());
            System.out.println();
            if (counter++ == 20 ) break;
            continue;
*/

            if ( ERR_MODE.equals("MM")) err = numMismatches(r,contig_seq);
            else if ( ERR_MODE.equals("MC") ) err = AlignmentUtils.numMismatches(r,contig_seq);
            else if ( ERR_MODE.equals("ERR")) err = numErrors(r,contig_seq);
            else if ( ERR_MODE.equals("MG")) err = numMismatchesGaps(r,contig_seq);
            if ( err > MAX_ERRS.intValue() ) continue;
            //        	counter++;
            //        	if ( counter % 1000000 == 0 ) System.out.println(counter+" records; "+col.memStatsString());
            col.receive(r);

        }
        
        if ( ! CONTROL_RUN ) {
            pileBuilder.printStats();
            pileBuilder.close();
        }
        System.out.println("done.");
        System.out.println("Discarded reads with non-M,I,D cigar elements: "+ discarded_cigar_count);
        System.out.println("Discarded long reads (above "+MAX_READ_LENGTH+" bp): "+ discarded_long_read_count);
        System.out.println();
        col.printLengthHistograms();
        samReader.close();
        ptWriter.close();
        return 0;
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

    /** Trivial utility method that goes some distance trying to ensure that the input file is there;
     * the only purpose is reducing clutter in main(). Receives a default
     * input file argument, does a few checks (e.g. that it is non-null and exists), if they fail tries
     * to fire up a file chooser dialog using start_folder as initial directory, etc.
     * @param default_arg some "default" input file; if it is non-null and exists, nothing else will be done,
     *        and the same default_arg objetc will be returned; otherwise the method will try to ask for a "better" input.
     * @param start_folder should file open dialog be fired up, it will initially display this directory.
     * @return File object that is not null and does exist (there is no check that it is a valid SAM/BAM file though).
     */
    private File getInputFile(File default_arg, String start_folder) {
        File f = default_arg;
        if ( f==null || ! f.exists() ) {
            JFileChooser fc = new JFileChooser(start_folder);
            FileNameExtensionFilter ff = new FileNameExtensionFilter("SAM and BAM files","sam","bam");
            fc.setFileFilter(ff);
            fc.setFileSelectionMode(JFileChooser.FILES_ONLY);

            int ret = fc.showOpenDialog(null);
            f = fc.getSelectedFile();
            if ( ret != JFileChooser.APPROVE_OPTION ) {
                System.out.println("No input file specified. Exiting...");
                System.exit(1);
            }
        }

        if ( f == null || ! f.exists() ) {
            System.out.println("SAM or BAM input file must be specified. Exiting...");
            System.exit(1);
        }

        return f;
    }

    /** Auxiliary method to remove some clutter from main(); gets called only once and tries to get
     * contig ordering from the header provided by opened SAM reader; if no header info is available
     * falls back to default ordering; whichever ordering is used, it is set for GenomeLoc class.
     * @param r sam reader to get header from
     */
    private void setContigOrdering(SAMFileReader r) {
        SAMFileHeader h = r.getFileHeader();
        if ( h == null ) {
            System.out.println("No header found in SAM file, falling back to default contig ordering");
            setDefaultContigOrdering();
            return;
        }
        GenomeLoc.setupRefContigOrdering(h.getSequenceDictionary());
    }

    private void setDefaultContigOrdering() {
        Map<String,Integer> rco = new HashMap<String,Integer>();
        rco.put("chrM",0);
        for ( int i = 1 ; i <= 22 ; i++ ) rco.put(Integer.toString(i),i);//rco.put("chr"+i,i);
        rco.put("chrX",23);
        rco.put("chrY",24);
    }
}
