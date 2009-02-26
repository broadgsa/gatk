/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

/**
 * Interface for SAMText and BAM file writers.  Clients need not care which they write to,
 * once the object is constructed.
 */
public interface SAMFileWriter {
    void addAlignment(SAMRecord alignment);

    /**
     * Must be called or file will likely be defective. 
     */
    void close();
}
