/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

/**
 * Optimized method for converting Solexa ASCII qualities into Phred scores.
 * Pre-computes all values in order to eliminate repeated computation.
 */
public class SolexaQualityConverter {

    /**
     * This value is added to a Solexa quality score to make it printable ASCII
     */
    private static int SOLEXA_ADDEND = 64;

    /**
     * Mapping from ASCII value in Gerald export file to phred score
     */
    private final byte[] phredScore = new byte[256];

    public SolexaQualityConverter() {
        for (int i = 0; i < SOLEXA_ADDEND; ++i) {
            phredScore[i] = 0;
        }
        for (int i = SOLEXA_ADDEND; i < phredScore.length; ++i) {
            phredScore[i] = decodeSolexaQualityToPhred(i);
        }
    }


    /** Converts a solexa character quality into a phred numeric quality. */
    private byte decodeSolexaQualityToPhred(final int solexaQuality) {
        return (byte) Math.round(10d * Math.log10(1d+Math.pow(10d, (solexaQuality - SOLEXA_ADDEND)/10d)));
    }

    /**
     * Convert a solexa quality ASCII character into a phred score.
     */
    public byte solexaToPhred(final byte solexaQuality) {
        return phredScore[solexaQuality];
    }

    /**
     * @return a byte array that can be indexed by Solexa ASCII quality, with value
     * of corresponding Phred score.  Elements 0-63 are invalid because Solexa qualities
     * should all be >= 64.  Do not modify this array!
     */
    public byte[] getSolexaToPhredConversionTable() {
        return phredScore;
    }
}
