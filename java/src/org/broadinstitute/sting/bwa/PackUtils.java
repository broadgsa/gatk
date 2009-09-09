package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.StingException;

/**
 * Utilities designed for packing / unpacking bases.
 *
 * @author mhanna
 * @version 0.1
 */
public class PackUtils {
    /**
     * How many possible bases can be encoded?
     */
    public static final int ALPHABET_SIZE = 4;

    /**
     * How many bits does it take to store a single base?
     */
    public static final int BITS_PER_BASE = (int)(Math.log(ALPHABET_SIZE)/Math.log(2));

    /**
     * How many bits fit into a single byte?
     */
    public static final int BITS_PER_BYTE = 8;

    /**
     * How many bits can a given type hold?
     * @param type Type to test.
     * @return Number of bits that the given type can hold.
     */
    public static final int bitsInType( Class<?> type ) {
        try {
            long typeSize = type.getField("MAX_VALUE").getLong(null) - type.getField("MIN_VALUE").getLong(null)+1;
            long intTypeSize = (long)Integer.MAX_VALUE - (long)Integer.MIN_VALUE + 1;
            if( typeSize > intTypeSize )
                throw new StingException("Cannot determine number of bits available in type: " + type.getName());
            return (int)(Math.log(typeSize)/Math.log(2));
        }
        catch( NoSuchFieldException ex ) {
            throw new StingException("Cannot determine number of bits available in type: " + type.getName(),ex);
        }
        catch( IllegalAccessException ex ) {
            throw new StingException("Cannot determine number of bits available in type: " + type.getName(),ex);
        }
    }

    /**
     * Gets the two-bit representation of a base.  A=00b, C=01b, G=10b, T=11b.
     * @param base ASCII value for the base to pack.
     * @return A byte from 0-3 indicating the base's packed value.
     */
    public static byte packBase(byte base) {
        switch( base ) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            default:
                throw new StingException("Unknown base type: " + base);
        }
    }

    public static byte unpackBase(byte pack) {
        switch( pack ) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                throw new StingException("Unknown pack type: " + pack);
        }
    }
}
