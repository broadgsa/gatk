package org.broadinstitute.sting.indels;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.HashMap;


import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequenceFileWalker;
import edu.mit.broad.picard.reference.ReferenceSequence;

import net.sf.samtools.*;
import org.broadinstitute.sting.utils.GenomeLoc;

public class IndelInspector extends CommandLineProgram {

    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Investigates indels called in the alignment data\n";
    @Option(shortName="I", doc="SAM or BAM file for calling") public File INPUT_FILE;
    @Option(shortName="L",doc="Genomic interval to run on, as contig[:start[-stop]]; whole genome if not specified", optional=true) public String GENOME_LOCATION;
    @Option(doc="Error counting mode: MM - count mismatches only, ERR - count errors (arachne style), MG - count mismatches and gaps as one error each") public String ERR_MODE;
    @Option(doc="Maximum number of errors allowed (see ERR_MODE)") public Integer MAX_ERRS;
//    @Option(shortName="R", doc="Reference fasta or fasta.gz file") public File REF_FILE;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new IndelInspector().instanceMain(argv));
    }
    
	protected int doWork() {

        GenomeLoc location = null;
        if ( GENOME_LOCATION != null ) {
            location = GenomeLoc.parseGenomeLoc(GENOME_LOCATION);
        }

		if ( ! ERR_MODE.equals("MM") && ! ERR_MODE.equals("MG") && ! ERR_MODE.equals("ERR") ) {
			System.out.println("Unknown value specified for ERR_MODE");
			return 1;
		}
		
		final SAMFileReader samReader = new SAMFileReader(getInputFile(INPUT_FILE,"/broad/1kG/"));

        setContigOrdering(samReader);

        ReferenceSequenceFileWalker reference = new ReferenceSequenceFileWalker(
                    new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
            );


        ReferenceSequence contig_seq = null;

		IndelRecordPileCollector col = null;
		try {
			col = new IndelRecordPileCollector();
		} catch(Exception e) { System.err.println(e.getMessage()); }
		if ( col == null ) return 1; 
		
        String cur_contig = null;
        int counter = 0;

        for ( SAMRecord r : samReader ) {

        	if ( r.getReferenceName() != cur_contig) {
        		cur_contig = r.getReferenceName();
        		System.out.println("Contig "+cur_contig);
                // if contig is specified and we are past that contig, we are done:
                if ( location != null && GenomeLoc.compareContigs(cur_contig, location.getContig()) == 1 ) break;
                if ( location == null || GenomeLoc.compareContigs(cur_contig, location.getContig()) == 0 ) {
                    contig_seq = reference.get(r.getReferenceIndex());
                    System.out.println("loaded contig "+cur_contig+" (index="+r.getReferenceIndex()+"); length="+contig_seq.getBases().length+" tst="+contig_seq.toString());
                }
            }

            // if contig is specified and wqe did not reach it yet, skip the records until we reach that contig:
        	if ( location != null && GenomeLoc.compareContigs(cur_contig, location.getContig()) == -1 ) continue;

            // if stop position is specified and we are past that, stop reading:
        	if ( location != null && r.getAlignmentStart() > location.getStop() ) break;
        	
        	if ( cur_contig.equals("chrM") || GenomeLoc.compareContigs(cur_contig,"chrY")==1 ) continue; // skip chrM and unplaced contigs for now
        	
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

        	if ( ERR_MODE.equals("MM")) err = numMismatches(r);
        	else if ( ERR_MODE.equals("ERR")) err = numErrors(r);
        	else if ( ERR_MODE.equals("MG")) err = numMismatchesGaps(r);
        	if ( err > MAX_ERRS.intValue() ) continue;
        	counter++;
        	if ( counter % 1000000 == 0 ) System.out.println(counter+" records; "+col.memStatsString());
        	col.receive(r);

        }
        System.out.println("done.");
        col.printLengthHistograms();
        samReader.close();
        return 0;
	}
	
	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD 
	 * 
	 * @param r SAM record that must specify an alignment
	 * @return number of errors (number of mismatches plus total length of all insertions/deletions
	 * @throws RuntimeException
	 */
	private static int numErrors(SAMRecord r) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
		int errs = numMismatches(r);
		
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
	private static int numMismatchesGaps(SAMRecord r) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
		int errs = numMismatches(r);
		
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

    private static int numMismatchesDirect(SAMRecord r, ReferenceSequence ref) {
        int i_ref = r.getAlignmentStart()-1; // position on the ref
        int i_read = 0;                    // position on the read
        int mm_count = 0; // number of mismatches
        Cigar c = r.getCigar();
        for ( int k = 0 ; k < c.numCigarElements() ; k++ ) {
            CigarElement ce = c.getCigarElement(k);
            switch( ce.getOperator() ) {
                case M:
                    for ( int l = 0 ; l < ce.getLength() ; l++ ) {
                        if ( Character.toUpperCase(r.getReadString().charAt(i_read) ) !=
                             Character.toUpperCase((char)ref.getBases()[i_ref]) ) mm_count++;
                        i_ref++;
                        i_read++;
                    }
                    break;
                case I: i_read += ce.getLength(); break;
                case D: i_ref += ce.getLength(); break;
                default: throw new RuntimeException("Unrecognized cigar element");
            }

        }
        return mm_count;
    }
	
	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD */
	private static int numMismatches(SAMRecord r) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
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
        List<SAMSequenceRecord> seqs = h.getSequences();
        if ( seqs == null ) {
            System.out.println("No reference sequence records found in SAM file header, " +
                               "falling back to default contig ordering");
            setDefaultContigOrdering();
            return;
        }
        int i = 0;
        Map<String,Integer> rco = new HashMap<String,Integer>();
        for ( SAMSequenceRecord sr : seqs) {
            rco.put(sr.getSequenceName(),i++);
        }
        GenomeLoc.setContigOrdering(rco);
    }

    private void setDefaultContigOrdering() {
        Map<String,Integer> rco = new HashMap<String,Integer>();
        rco.put("chrM",0);
        for ( int i = 1 ; i <= 22 ; i++ ) rco.put("chr"+i,i);
        rco.put("chrX",23);
        rco.put("chrY",24);
    }
}
