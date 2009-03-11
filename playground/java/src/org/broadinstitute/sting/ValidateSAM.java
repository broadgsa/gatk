package org.broadinstitute.sting;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.cmdline.Option;

import java.io.*;

public class ValidateSAM extends CommandLineProgram {
    // Usage and parameters
    @Usage(programVersion="0.1") public String USAGE = "SAM Validator\n";
    @Option(shortName="I", doc="SAM or BAM file for validation") public File INPUT_FILE;
    @Option(shortName="M", doc="Maximum number of errors to detect before exiting", optional=true) public String MAX_ERRORS_ARG = "-1";
    @Option(shortName="S", doc="How strict should we be with validation", optional=true) public String STRICTNESS_ARG = "strict";

    private long startTime = -1;
    
    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new ValidateSAM().instanceMain(argv));
    }

    public void printProgress( int nRecords, int nErrors ) {
        final double elapsed = (System.currentTimeMillis() - startTime) / 1000.0;
        final double secsPer1MReads = (elapsed * 1000000.0) / nRecords;
    	System.out.printf("Read %d records containing %d errors in %.2f secs (%.2f secs per 1M reads)%n", nRecords, nErrors, elapsed, secsPer1MReads);
    }
    
    protected int doWork() {
    	int MAX_ERRORS = -1;	// Don't bail ever
        if ( MAX_ERRORS_ARG != null ) {
        	MAX_ERRORS = Integer.parseInt(MAX_ERRORS_ARG);
    	}
    	
        // Start the timer
    	startTime = System.currentTimeMillis();

    	// Initialize the sam reader
    	CloseableIterator<SAMRecord> iter = null;
      	try {
      		 final SAMFileReader samReader = getSamReader(INPUT_FILE);
      		 iter = samReader.iterator();
        } catch (Exception ioe) {
        	System.out.println("[VALIDATION FAILURE IN HEADER]: " + ioe);
            ioe.printStackTrace();
            return 1;
        }
        
        int nRecords = 0;
        int nErrors = 0;
        while ( iter.hasNext() ) {
        	nRecords++;
        	try {
            	final SAMRecord ri = iter.next();
	        } catch (Exception ioe) {
	        	nErrors++;
	        	System.out.println("[VALIDATION FAILURE IN RECORD]: " + ioe);
                ioe.printStackTrace();
	        }
	        
	        if ( MAX_ERRORS > -1 && nErrors >= MAX_ERRORS ) {
	        	System.out.println("Maximum number of errors encountered " + nErrors);
	        	break;
	        }
	        
	        if ( nRecords % 100000 == 0 ) {
	        	printProgress( nRecords, nErrors );
	        }
        }

    	printProgress( nRecords, nErrors );
        return 0;
    }

    private static void usage() {
        System.err.println("USAGE: org.broadinstitute.sting.ValidateSAM <SAMFile|BAMFile>");
    }

    private SAMFileReader getSamReader(final File samFile) {
   	
    	ValidationStringency strictness = SAMFileReader.ValidationStringency.STRICT;
    	if ( STRICTNESS_ARG == null ) {
            strictness = SAMFileReader.ValidationStringency.STRICT;
    	}
    	else if ( STRICTNESS_ARG.toLowerCase().equals("lenient") ) {
    		strictness = SAMFileReader.ValidationStringency.LENIENT; 
    	}
    	else if ( STRICTNESS_ARG.toLowerCase().equals("silent") ) {
    		strictness = SAMFileReader.ValidationStringency.SILENT;
    	}
    	else {
            strictness = SAMFileReader.ValidationStringency.STRICT;
    	}
            
        System.err.println("Strictness is " + strictness);   	
    	final SAMFileReader samReader = new SAMFileReader(samFile, true);
    	samReader.setValidationStringency(strictness);
        
        return samReader;
    }

}
