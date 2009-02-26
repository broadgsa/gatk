/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.sam;

/**
 * Constants for tags used in our SAM/BAM files
 */
public class ReservedTagConstants {
    public static final String READ_GROUP_ID = "RG"; // Specified in the SAM spec doc
    public static final String XN = "XN";    // Present and set to 1 if a read is a noise read
}
