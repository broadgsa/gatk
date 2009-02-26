/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.aligner.maq;

/**
 * Utility class to hold Maq-related constants (program name, location, switches, etc)
 */
public class MaqConstants {
    // General Maq constants
    public static final String PROGRAM_NAME = "Maq";
    public static final String PROGRAM_VERSION = "0.7.1";
    public static final String MAQ_HOME = "/seq/dirseq/maq-0.7.1/";

    // Command-related constants
    public static final String MAQ_COMMAND = "maq";
    public static final String MAP_COMMAND = "map";
    public static final String MERGE_COMMAND = "mapmerge";

    // Constants related to Maq map switches
    public static final String SWITCH_SUM_MISMATCHES = "-e";
    public static final int HIGH_STRINGENCY_SUM_MISMATCHES = 100;
    public static final int LOW_STRINGENCY_QUALITY_FOR_MISMATCHES = 30;

    public static final String SWITCH_MAX_OUTER_DISTANCE = "-a";
    public static final int LOW_STRINGENCY_MAX_OUTER_DISTANCE = 1500;
    public static final double HIGH_STRINGENCY_MAX_OUTER_DISTANCE_MULTIPLIER = 1.5d;

    public static final String SWITCH_RANDOM_SEED =   "-s";
    public static final int DEFAULT_RANDOM_SEED = 0;

    public static String getProgramVersion() { return PROGRAM_VERSION; }
}
