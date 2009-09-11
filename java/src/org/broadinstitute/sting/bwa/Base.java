package org.broadinstitute.sting.bwa;

import java.util.EnumSet;
import java.util.Map;
import java.util.HashMap;

/**
 * Enhanced enum representation of a base.
 *
 * @author mhanna
 * @version 0.1
 */
public enum Base
{
    A((byte)'A',0),
    C((byte)'C',1),
    G((byte)'G',2),
    T((byte)'T',3);

    /**
     * The ASCII representation of a given base.
     */
    private final byte ascii;

    /**
     * The 2-bit packed value of the base.
     */
    private final int pack;

    /**
     * Representation of the base broken down by packed value.
     */
    private static final Map<Integer,Base> basesByPack = new HashMap<Integer,Base>();

    /**
     * Representation of the base broken down by ASCII code.
     */
    private static final Map<Byte,Base> basesByASCII = new HashMap<Byte,Base>();

    static {
        for(Base base : EnumSet.allOf(Base.class)) {
            basesByPack.put(base.pack,base);
            basesByASCII.put(base.ascii,base);
        }
    }

    /**
     * Create a new base with the given ascii representation and
     * pack value.
     * @param ascii ASCII representation of a given base.
     * @param pack Packed value of a given base.
     */
    private Base( byte ascii, int pack ) {
        this.ascii = ascii;
        this.pack = pack;
    }

    /**
     * Get the given base from the packed representation.
     * @param pack Packed representation.
     * @return base.
     */
    public static Base fromPack( int pack ) { return basesByPack.get(pack); }

    /**
     * Convert the given base to its packed value.
     * @return Packed value.
     */
    public int toPack() { return pack; }

    /**
     * Convert the given base to its packed value.
     * @param ascii ASCII representation of the base.
     * @return Packed value.
     */
    public static int toPack( byte ascii ) { return basesByASCII.get(ascii).pack; }

    /**
     * Get the given base from the ASCII representation.
     * @param ascii ASCII representation.
     * @return base.
     */
    public static Base fromASCII( byte ascii ) { return basesByASCII.get(ascii); }    

    /**
     * Convert the given base to its ASCII value.
     * @return ASCII value.
     */
    public byte toASCII() { return ascii; }

    /**
     * Convert the given base to its ASCII value.
     * @param pack The packed representation of the base.
     * @return ASCII value.
     */
    public static byte toASCII( int pack ) { return basesByPack.get(pack).ascii; }
}
