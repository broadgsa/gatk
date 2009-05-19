package org.broadinstitute.sting.playground.indels;

import java.io.File;

import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequenceFileWalker;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class ShowMSA extends CommandLineProgram {

    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Prints MSA into stdout\n";
    @Option(shortName="I", doc="SAM or BAM file with alignment data") public File INPUT_FILE;
    @Option(shortName="L", doc="Contig:Start-Stop or Contig:poslocation of the window to draw") public String LOCATION;
    @Option(shortName="W", doc="Number of bases on each side of specified position if LOCATION is in Contig:pos format; ignored otherwise", optional=true) public Integer WINDOW;
    @Option(shortName="R", doc="Reference fastb file") public File REF_FILE;
    @Option(shortName="P", doc="If true, then any read (partially) overlapping with the specified region will be shown. "+
    		"Otherwise (default), only reads fully contained in the specified interval are shown", optional=true) public Boolean PARTIAL;
    @Option(doc="Error counting mode: MM - count mismatches only, ERR - count errors (arachne style), MG - count mismatches and gaps as one error each") public String ERR_MODE;
    @Option(doc="Maximum number of errors allowed (see ERR_MODE)") public Integer MAX_ERRS;
    @Option(shortName="F",doc="Format: PILE - show alignment, FASTA - print sequences in fasta",optional=true) public String OUT_FORMAT;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new ShowMSA().instanceMain(argv));
    }
    
	protected int doWork() {
		
		if ( ! ERR_MODE.equals("MM") && ! ERR_MODE.equals("MG") && ! ERR_MODE.equals("ERR") ) {
			System.out.println("Unknown value specified for ERR_MODE");
			return 1;
		}

		if ( PARTIAL == null ) PARTIAL = new Boolean(false);
		if ( OUT_FORMAT == null ) OUT_FORMAT=new String("PILE");
		
		if ( ! OUT_FORMAT.equals("PILE") && ! OUT_FORMAT.equals("FASTA")) {
			System.out.println("OUT_FORMAT can only have values PILE or FASTA");
			return 1;
		}
		
		if ( ! INPUT_FILE.exists() ) {
			System.out.println("Specified INPUT_FILE does not exist");
			return 1;
		}

		if ( ! REF_FILE.exists() ) {
			System.out.println("Specified REF_FILE does not exist");
			return 1;
		}

		if ( LOCATION.indexOf(':') == -1 ) {
			System.out.println("LOCATION should follow Contig:Start-Stop or Contig:Pos format");
			return 1;
		}
		String[] s1 = LOCATION.split(":");
		int contig;
		try {
			contig = Integer.valueOf(s1[0]);
		}	catch (NumberFormatException e) {
			System.out.println("LOCATION: contig must be specified as an integer");
			return 1;
		}
		
		if ( s1.length != 2 ) {
			System.out.println("LOCATION should follow Contig:Start-Stop or Contig:Pos format");
			return 1;
		}
		
		String s2[] = s1[1].split("-");
		if ( s2.length > 2 ) {
			System.out.println("LOCATION should follow Contig:Start-Stop or Contig:Pos format");
			return 1;
		}
		int left, right;
		if ( s2.length == 2 ) {
			try {
				left = Integer.valueOf(s2[0]);
				right = Integer.valueOf(s2[1]);
			}	catch (NumberFormatException e) {
				System.out.println("LOCATION: window boundaries should be specified as integers");
				return 1;
			}
		} else {
			int pos = 0;
			try {
				pos = Integer.valueOf(s2[0]);
			}	catch (NumberFormatException e) {
				System.out.println("LOCATION: position on the contig should be specified as an integer");
				return 1;
			}
			if (WINDOW == null ) {
				System.out.println("WINDOW must be specified when LOCATION specifies a single poisiton (Contig:Pos)");
				return 1;
			}
			left = pos - WINDOW.intValue();
			right = pos+WINDOW.intValue();
		}
		
		
		String ref_contig ;
		
		try {
		    ReferenceSequenceFileWalker mRefReader =
                    new ReferenceSequenceFileWalker(ReferenceSequenceFileFactory.getReferenceSequenceFile(REF_FILE));
			ref_contig = mRefReader.get(contig).toString(); // reload ref
		} catch (Exception e) {
			System.out.println("Failed to read reference sequence from " + REF_FILE);
			return 1;
		}

		SAMFileReader reader ;
		try {
			reader = new SAMFileReader(INPUT_FILE);
		} catch ( Exception e) {
			System.out.println(e.getMessage());
			return 1;
		}
		
		SequencePile msa=null;
		
		if ( OUT_FORMAT.equals("PILE")) {
			msa = new SequencePile(ref_contig.substring(left-1, right));
		} else {
			System.out.println(">reference "+contig+":"+left+"-"+right);
			System.out.println(ref_contig.substring(left-1, right));
		}
		
		for( SAMRecord r : reader ) {
			if ( r.getReadUnmappedFlag() ) continue;
			if ( r.getReferenceIndex() < contig ) continue;
			if ( r.getReferenceIndex() > contig ) break;
			if ( r.getAlignmentEnd() < left ) continue;
			if ( r.getAlignmentStart() >= right ) break;
			if ( ! PARTIAL && ( r.getAlignmentStart() < left || r.getAlignmentEnd() >= right ) ) continue;
			
        	int err = -1;
        	if ( ERR_MODE.equals("MM")) err = numMismatches(r);
        	else if ( ERR_MODE.equals("ERR")) err = numErrors(r);
        	else if ( ERR_MODE.equals("MG")) err = numMismatchesGaps(r);
        	if ( err > MAX_ERRS ) continue;

        	if ( OUT_FORMAT.equals("PILE") ) {
        		msa.addAlignedSequence(r.getReadString(), r.getReadNegativeStrandFlag(), r.getCigar(), r.getAlignmentStart() - left);
        	} else {
        		System.out.print(">read "+r.getReadName());
        		if ( r.getReadNegativeStrandFlag() ) System.out.println("(rc)");
        		else System.out.println("(fw)");
        		System.out.println(r.getReadString());
        	}
		}
		
		if ( OUT_FORMAT.equals("PILE") ) msa.colorprint();
////			System.out.println(msa.format());
		
		return 0;
	}

	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD 
	 * 
	 * @param r SAM record that must specify an alignment
	 * @return number of errors (number of mismatches plus total length of all insertions/deletions
	 * @throws RuntimeException if cigar contains any elements other than M,I,D
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
	 * @throws RuntimeException if cigar contains any elements other than M,I,D
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
	
	
	/** This method is a HACK: it is designed to work around the current bug in NM tags created  at CRD */
	private static int numMismatches(SAMRecord r) throws RuntimeException {
		
		// NM currently stores the total number of mismatches in all blocks + 1
		return ((Integer)r.getAttribute("NM")).intValue() - 1;
		
	}
	
}
