package org.broadinstitute.sting.utils.clipreads;

/**
 * How should we represent a clipped bases in a read?
 */
public enum ClippingRepresentation {
    WRITE_NS,           // change the bases to Ns
    WRITE_Q0S,          // change the quality scores to Q0
    WRITE_NS_Q0S,       // change the quality scores to Q0 and write Ns
    SOFTCLIP_BASES,     // change cigar string to S, but keep bases
    HARDCLIP_BASES      // remove the bases from the read
}
