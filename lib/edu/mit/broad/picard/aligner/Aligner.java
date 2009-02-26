/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.aligner;

/**
 * API for aligners.  Clients must call these methods in order, as each depends on
 * the previous one, but they may call them multiple times and need not call them all.
 * This allows steps to be rerun and also lets the caller review intermediate files 
 * when troubleshooting.
 *
 * @author Kathleen Tibbetts
 */
public interface Aligner {

    public static enum Stringency{ low, high };

    /**
     * Prepares all the necessary inputs for the alignment process from a BAM file of read data.
     */
    public void prepareInputs();

    /**
     * Does the alignment and produces output in the underlying form of the aligner.
     */
    public void align();

    /**
     * Converts the output of the aligner to BAM format
     */
    public void prepareOutput();

    /**
     * Cleans up intermediate files (the files created in by and for the underlying aligner by the
     * prepareInputs() and align() methods.  Does not clean up the original source files or the final BAM file.
     */
    public void cleanup();

}
