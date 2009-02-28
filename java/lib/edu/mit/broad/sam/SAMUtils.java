/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.sam;


/**
 * Utilty methods.
 */
final class SAMUtils
{
    private static final byte COMPRESSED_EQUAL_LOW = 0;
    private static final byte COMPRESSED_A_LOW = 1;
    private static final byte COMPRESSED_C_LOW = 2;
    private static final byte COMPRESSED_G_LOW = 4;
    private static final byte COMPRESSED_T_LOW = 8;
    private static final byte COMPRESSED_N_LOW = 15;
    private static final byte COMPRESSED_EQUAL_HIGH = COMPRESSED_EQUAL_LOW << 4;
    private static final byte COMPRESSED_A_HIGH = COMPRESSED_A_LOW << 4;
    private static final byte COMPRESSED_C_HIGH = COMPRESSED_C_LOW << 4;
    private static final byte COMPRESSED_G_HIGH = COMPRESSED_G_LOW << 4;
    private static final byte COMPRESSED_T_HIGH = (byte)(COMPRESSED_T_LOW << 4);
    private static final byte COMPRESSED_N_HIGH = (byte)(COMPRESSED_N_LOW << 4);

    private SAMUtils() {
    }

    static int unpackInt16(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8));
    }

    static int unpackInt32(final byte[] buffer, final int offset) {
        return ((buffer[offset] & 0xFF) |
                ((buffer[offset+1] & 0xFF) << 8) |
                ((buffer[offset+2] & 0xFF) << 16) |
                ((buffer[offset+3] & 0xFF) << 24));
    }

    /**
     * Convert from a byte array containing =AaCcGgTtNn, to a byte array half as long,
     * with =, A, C, G, T converted to 0, 1, 2, 4, 8, 15
     * @param readBases
     * @return
     */
    static byte[] bytesToCompressedBases(final byte[] readBases) {
        final byte[] compressedBases = new byte[(readBases.length + 1)/2];
        int i;
        for (i = 1; i < readBases.length; i+=2) {
            compressedBases[i/2] = (byte)(charToCompressedBaseHigh(readBases[i-1]) |
                                    charToCompressedBaseLow(readBases[i]));
        }
        // Last nybble
        if (i == readBases.length) {
            compressedBases[i/2] = charToCompressedBaseHigh((char)readBases[i-1]);
        }
        return compressedBases;
    }

    static byte[] compressedBasesToBytes(final int length, final byte[] compressedBases, final int compressedOffset) {
        final byte[] ret = new byte[length];
        int i;
        for (i = 1; i < length; i+=2) {
            ret[i-1] = compressedBaseToByteHigh(compressedBases[i/2 + compressedOffset]);
            ret[i] = compressedBaseToByteLow(compressedBases[i/2 + compressedOffset]);
        }
        // Last nybble
        if (i == length) {
            ret[i-1] = compressedBaseToByteHigh(compressedBases[i/2 + compressedOffset]);
        }
        return ret;
    }

    /**
     *
     * @param base One of =AaCcGgTtNn
     * @return nybble-encoded equivalent
     */
    private static byte charToCompressedBaseLow(final int base) {
        switch (base) {
            case '=':
                return COMPRESSED_EQUAL_LOW;
            case 'a':
            case 'A':
                return COMPRESSED_A_LOW;
            case 'c':
            case 'C':
                return COMPRESSED_C_LOW;
            case 'g':
            case 'G':
                return COMPRESSED_G_LOW;
            case 't':
            case 'T':
                return COMPRESSED_T_LOW;
            case 'n':
            case 'N':
            case '.':
                return COMPRESSED_N_LOW;
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    private static byte charToCompressedBaseHigh(final int base) {
        switch (base) {
            case '=':
                return COMPRESSED_EQUAL_HIGH;
            case 'a':
            case 'A':
                return COMPRESSED_A_HIGH;
            case 'c':
            case 'C':
                return COMPRESSED_C_HIGH;
            case 'g':
            case 'G':
                return COMPRESSED_G_HIGH;
            case 't':
            case 'T':
                return COMPRESSED_T_HIGH;
            case 'n':
            case 'N':
            case '.':
                return COMPRESSED_N_HIGH;
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    /**
     *
     * @param base One of COMPRESSED_*
     * @return one of ACGTN=
     */
    private static byte compressedBaseToByteLow(final int base) {
        switch (base & 0xf) {
            case COMPRESSED_EQUAL_LOW:
                return '=';
            case COMPRESSED_A_LOW:
                return 'A';
            case COMPRESSED_C_LOW:
                return 'C';
            case COMPRESSED_G_LOW:
                return 'G';
            case COMPRESSED_T_LOW:
                return 'T';
            case COMPRESSED_N_LOW:
                return 'N';
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    private static byte compressedBaseToByteHigh(final int base) {
        switch ((byte)(base & 0xf0)) {
            case COMPRESSED_EQUAL_HIGH:
                return '=';
            case COMPRESSED_A_HIGH:
                return 'A';
            case COMPRESSED_C_HIGH:
                return 'C';
            case COMPRESSED_G_HIGH:
                return 'G';
            case COMPRESSED_T_HIGH:
                return 'T';
            case COMPRESSED_N_HIGH:
                return 'N';
            default:
                throw new IllegalArgumentException("Bad  byte passed to charToCompressedBase: " + base);
        }
    }

    static String bytesToHexString(final byte[] data) {
        final char[] chars = new char[2 * data.length];
        for (int i = 0; i < data.length; i++) {
            final byte b = data[i];
            chars[2*i] = toHexDigit((b >> 4) & 0xF);
            chars[2*i+1] = toHexDigit(b & 0xF);
        }
        return new String(chars);
    }

    static byte[] hexStringToBytes(final String s)  throws NumberFormatException {
        if (s.length() % 2 != 0) {
            throw new NumberFormatException("Hex representation of byte string does not have even number of hex chars: " + s);
        }
        final byte[] ret = new byte[s.length() / 2];
        for (int i = 0; i < ret.length; ++i) {
            ret[i] = (byte) (fromHexDigit(s.charAt(i * 2)) << 4 + fromHexDigit(s.charAt(i * 2 + 1)));
        }
        return ret;
    }

    static String phredToFastq(final byte[] data) {
        if (data == null) {
            return null;
        }
        return phredToFastq(data, 0, data.length);
    }

    static String phredToFastq(final byte[] buffer, final int offset, final int length) {
        final char[] chars = new char[length];
        for (int i = 0; i < length; i++) {
            chars[i] = phredToFastq(buffer[offset+i] & 0xFF);
        }
        return new String(chars);
    }

    static char phredToFastq(final int phredScore) {
        if (phredScore < 0 || phredScore > 63) {
            throw new IllegalArgumentException("Cannot encode phred score: " + phredScore);
        }
        return (char) (33 + phredScore);
    }

    static byte[] fastqToPhred(final String fastq) {
        if (fastq == null) {
            return null;
        }
        final int length = fastq.length();
        final byte[] scores = new byte[length];
        for (int i = 0; i < length; i++) {
            scores[i] = (byte) fastqToPhred(fastq.charAt(i));
        }
        return scores;
    }

    static int fastqToPhred(final char ch) {
        if (ch < 33 || ch > 126) {
            throw new IllegalArgumentException("Invalid fastq character: " + ch);
        }
        return (ch - 33);
    }

    private static char toHexDigit(final int value) {
        return (char) ((value < 10) ? ('0' + value) : ('A' + value - 10));
    }

    private static int fromHexDigit(final char c) throws NumberFormatException {
        final int ret = Character.digit(c, 16);
        if (ret == -1) {
            throw new NumberFormatException("Not a valid hex digit: " + c);
        }
        return ret;
    }

    /**
     * calculate the bin given an alignment in [beg,end)
     * Copied from SAM spec. 
     */
    static int reg2bin(final int beg, int end)
    {

        --end;

        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return  ((1<<9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return  ((1<<6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return  ((1<<3)-1)/7 + (beg>>26);
        return 0;
    }
}
