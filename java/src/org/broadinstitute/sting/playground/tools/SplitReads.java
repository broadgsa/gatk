package org.broadinstitute.sting.playground.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.io.File;

import net.sf.samtools.*;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;

public class SplitReads extends CommandLineProgram {
    @Usage(programVersion="1.0") public String USAGE = "Splits reads: extracts sub-sequences of the specified length(s) from left "+
            "and/or right ends of all the reads into the specified output bam file(s). For the reads in the input that are mapped, "+
            "the subsequences in the output bam(s) will have appropriately adjusted alignment positions and chopped cigars.";
    @Option(shortName="I",
            doc="Input file (bam or sam) with read sequences to split.",
            optional=false)
    public File IN = null;
    @Option(shortName="E", doc="Read end to select, 1=left, 2=right; default: select both ends.",
    		optional=true) public List<Integer> READ_ENDS = new ArrayList<Integer>();
    @Option(shortName="N", doc="Number of bases to keep in the corresponding segment of the read. "+
    			"Synchronized with READ_ENDS argument; if single number is given, all selected segments (ends) will have specified length.",
    		optional=false) public List<Integer> LENGTH = new ArrayList<Integer>();
    @Option(shortName="S", doc="Read name for each segment (read end) will be set as original read name followed by the corresponding suffix." +
    			"Synchronized with READ_ENDS argument and must have the same number of entries if specified (note that default READ_ENDS is a list of (1,2). "+
    			"By default, suffixes are empty strings, i.e. all segments have the same name(s) as the original read." , optional=true) public List<String> SUFFIXES = new ArrayList<String>();
    @Option(shortName="O",optional=false, doc="Each read end will be sent into the corresponding file " +
    		"(synchronized with READ_ENDS). If only one file name is specified, all read segments will be printed into that file."
    		) public List<File> OUTPUT_BAMS = new ArrayList<File>();
    @Option(shortName="U", doc="Split and output only unmapped reads; mapped reads will be ignored.",
    		optional=true) public boolean UNMAPPED = false;


    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new SplitReads().instanceMain(argv));
    }

    protected int doWork() {

    	// if read ends are not specified explicitly on the cmd line, set default 1,2 (both ends)
    	if ( READ_ENDS.size() == 0 ) {
    		READ_ENDS.add(1);
    		READ_ENDS.add(2);
    	}

    	for ( Integer i : READ_ENDS) {
    		if ( ! i.equals(1) && ! i.equals(2)) throw new RuntimeException("Unknown value specified for READ_ENDS: "+i);
    	}
    	
    	// if suffixes are not specified, set them to "", ""
    	if ( SUFFIXES.size() == 0 ) {
    		for ( Integer i : READ_ENDS) {
    			SUFFIXES.add( "" );
    		}
    	} else {
    		// or make sure that the number of suffixes matches the number of ends
    		if ( SUFFIXES.size() != READ_ENDS.size() ) throw new RuntimeException("Number of suffixes specified must be equal to the number of read ends requested."+
    				"Passed: "+ READ_ENDS.size() +" READ_ENDS and " + SUFFIXES.size() + " SUFFIXES arguments.");
    	}
    	
    	if ( LENGTH.size() == 1 ) {
    		// if only one length is specified, apply it to all ends:
    		LENGTH = Collections.nCopies(READ_ENDS.size(), LENGTH.get(0));
    	}
    	
    	if ( LENGTH.size() != READ_ENDS.size() ) throw new RuntimeException("Number of lengths specified must be equal to the number of read ends requested."+
				"Passed: "+ READ_ENDS.size() +" READ_ENDS and " + LENGTH.size() + " LENGTH arguments.");
    	
    	if ( READ_ENDS.size() != OUTPUT_BAMS.size() && OUTPUT_BAMS.size() != 1 )
            throw new RuntimeException("Number of output files must be either one, or equal to the number of read ends requested."+
				"Passed: "+ READ_ENDS.size() +" READ_ENDS and " + OUTPUT_BAMS.size() + " OUTPUT_BAMS arguments.");

        SAMFileReader inReader = new SAMFileReader(IN);

        List<SAMFileWriter> outWriters = new ArrayList<SAMFileWriter>(OUTPUT_BAMS.size());
        for ( File outName : OUTPUT_BAMS ) {
            outWriters.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(inReader.getFileHeader(), true, outName)) ;
        }


        for ( SAMRecord read : inReader ) {

            if ( UNMAPPED && ! read.getReadUnmappedFlag() ) continue;
            
            for ( int i = 0 ; i < READ_ENDS.size(); i++ ) {

                SAMRecord newRecord = null;
                try {
                    newRecord = (SAMRecord)read.clone();
                } catch (CloneNotSupportedException e) {
                    throw new RuntimeException("Clone not supported by SAMRecord implementation");
                }

                final int whichEnd = READ_ENDS.get(i);
                final int length = LENGTH.get(i);
                String name = read.getReadName();
                if ( length > read.getReadLength() ) throw new RuntimeException("Read "+name+" is shorter than the specified length ("+read.getReadLength()+"<"+length+")");
                int start = 0 , stop = 0; // [start, stop) : segment of the read to be selected; coordinates are wrt read sequence; half-open 0 based
                switch ( whichEnd ) {
                case 1: start = 0 ; stop = start + LENGTH.get(i); break;
                case 2: stop = read.getReadLength() ; start = stop - LENGTH.get(i); break;
                }

                newRecord.setReadBases(Arrays.copyOfRange(read.getReadBases(),start,stop));
                newRecord.setBaseQualities(Arrays.copyOfRange(read.getBaseQualities(), start, stop));
                newRecord.setReadName(name+ SUFFIXES.get(i));
                if ( read.getReadUnmappedFlag() ) {
                    //newRecord.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                } else {
                    newRecord.setAlignmentStart(read.getAlignmentStart()+start);
                    newRecord.setCigar( chopCigar(read.getCigar(), start, length ));
                }

                if ( outWriters.size() > 1 ) outWriters.get(i).addAlignment(newRecord);
                else outWriters.get(0).addAlignment(newRecord);
            }

        }

        inReader.close();
        for ( SAMFileWriter w : outWriters ) w.close();

        return 0;
    }
    

	/**
	 * Returns new cigar representing segment of the alignment that starts at position <code>start</code> (0-based)
	 * with respect to the start of the original cigar and covers <code>length</code> bases on the original read the
	 * <code>origCigar</code> corresponds to (i.e. I elements count, but D do not).
	 * @param origCigar
	 * @param start
	 * @param length
	 * @return
	 */
	private Cigar chopCigar( Cigar origCigar, int start, int length ) {

		int elementEnd = 0; // next base after the end of the current cigar element on the read
		
		Cigar newCigar = new Cigar();
		
		Iterator<CigarElement> elements = origCigar.getCigarElements().iterator();

        if ( ! elements.hasNext() ) System.out.println("CIGAR HAS NO ELEMENTS!");

		CigarElement ce = null;
		
		while ( elementEnd <= start ) { // if we did not reach the start of selected segment yet:
//            System.out.println("INIT: start="+start+"; length="+length+"; elementEnd="+elementEnd);
			ce = elements.next();
			switch ( ce.getOperator() ) {
			case N:  //
			case D : // read misses bases wrt the ref, nothing to count on the read
				break;
			case I:
			case M:
			case S:
			case H: // all these elements are real bases on the read. Skip them completely if 
				    // 'start' is past them, or crop if it is inside:
				elementEnd += ce.getLength(); // 1 base past end of the current element on the read

			}
		}
		// at this point we are guaranteed that ce is the element that contains 'start' position;
		// now we start adding cigar elements:

		// add manually first element, since we need only a part of it after 'start':
		newCigar.add( new CigarElement(Math.min(elementEnd-start, length), ce.getOperator()) );

		int selectionEnd = start + length;
//		    System.out.println(origCigar.toString()+": start="+start+"; length="+length+"; selectionEnd="+selectionEnd+"; elementEnd="+elementEnd);
		while ( elementEnd < selectionEnd ) {
			ce = elements.next();
			switch ( ce.getOperator() ) {
			case N:  //
			case D : // read misses bases wrt the ref, nothing to count on the read, but the element has to be added:
				newCigar.add( new CigarElement(ce.getLength(), ce.getOperator()) );				
				break;
			case I:
			case M:
			case S:
			case H: // all these elements are real bases on the read. Add them and count them 
				    // making sure that the last element gets cropped if needed:
				elementEnd += ce.getLength(); // 1 base past end of the current element on the read
				if ( elementEnd > selectionEnd ) { // this is the last element we have to consider and it needs to be cropped:
					newCigar.add( new CigarElement(ce.getLength() -  elementEnd + selectionEnd , ce.getOperator()) );									
				} else {
					newCigar.add( new CigarElement(ce.getLength(), ce.getOperator()) );									
				}
			}
			
		}
		return newCigar;
		
	}
	
}
