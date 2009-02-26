/*
  The Broad Institute
  SOFTWARE COPYRIGHT NOTICE AGREEMENT
  This software and its documentation are copyright 2009 by the
  Broad Institute/Massachusetts Institute of Technology. All rights are
  reserved.

  This software is supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
  use, misuse, or functionality.
*/
package edu.mit.broad.sam;

/**
 * Convert between string and internal CIGAR representations
 */
public class TextCigarCodec
{
    private static final byte ZERO_BYTE = "0".getBytes()[0];
    private static final byte NINE_BYTE = "9".getBytes()[0];
    
    private static final TextCigarCodec singleton = new TextCigarCodec();

    /**
     * It is not necssary to get the singleton but it is preferrable to use the same one
     * over and over vs. creating a new object for each BAMRecord.
     */
    static TextCigarCodec getSingleton() {
        return singleton;
    }


    /**
     * Convert from interal CIGAR representation to String
     */
    String encode(final Cigar cigar) {
        if (cigar.numCigarElements() == 0) {
            return SAMRecord.NO_ALIGNMENT_CIGAR;
        }
        final StringBuilder ret = new StringBuilder();
        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            ret.append(cigarElement.getLength());
            ret.append(cigarElement.getOperator());
        }
        return ret.toString();
    }

    Cigar decode(final String textCigar) {
        if (SAMRecord.NO_ALIGNMENT_CIGAR.equals(textCigar)) {
            return new Cigar();
        }
        final Cigar ret = new Cigar();
        final byte[] cigarBytes = textCigar.getBytes();
        for (int i = 0; i < cigarBytes.length; ++i) {
            if (!isDigit(cigarBytes[i])) {
                throw new IllegalArgumentException("Malformed CIGAR string: " + textCigar);
            }
            int length = (cigarBytes[i] - ZERO_BYTE);
            for (++i; isDigit(cigarBytes[i]); ++i) {
                length = (length * 10) + cigarBytes[i] - ZERO_BYTE;
            }
            final CigarOperator operator = CigarOperator.characterToEnum(cigarBytes[i]);
            ret.add(new CigarElement(length, operator));
        }
        return ret;
    }
    
    private boolean isDigit(final byte c) {
        return c >= ZERO_BYTE && c <= NINE_BYTE;
    }

    
        
}

/******************************************************************/
/**************************[END OF TextCigarCodec.java]*************************/
/******************************************************************/
