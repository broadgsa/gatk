package org.broadinstitute.sting.alignment.reference.bwt;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Enhanced enum representation of a base.
 *
 * @author mhanna
 * @version 0.1
 */
public class Bases implements Iterable<Byte>
{
    public static byte A = 'A';
    public static byte C = 'C';
    public static byte G = 'G';
    public static byte T = 'T';

    public static final Bases instance = new Bases();

    private static final List<Byte> allBases;

    /**
     * Representation of the base broken down by packed value.
     */
    private static final Map<Integer,Byte> basesByPack = new HashMap<Integer,Byte>();

    static {
        List<Byte> bases = new ArrayList<Byte>();
        bases.add(A);
        bases.add(C);
        bases.add(G);
        bases.add(T);
        allBases = Collections.unmodifiableList(bases);

        for(int i = 0; i < allBases.size(); i++)
            basesByPack.put(i,allBases.get(i));
    }

    /**
     * Create a new base with the given ascii representation and
     * pack value.
     */
    private Bases() {
    }

    /**
     * Return all possible bases.
     * @return Byte representation of all bases.
     */
    public static Collection<Byte> allOf() {
        return allBases;
    }

    /**
     * Gets the number of known bases.
     * @return The number of known bases.
     */
    public static int size() {
        return allBases.size();
    }

    /**
     * Gets an iterator over the total number of known base types.
     * @return Iterator over all known bases.
     */
    public Iterator<Byte> iterator() {
        return basesByPack.values().iterator();
    }

    /**
     * Get the given base from the packed representation.
     * @param pack Packed representation.
     * @return base.
     */
    public static byte fromPack( int pack ) { return basesByPack.get(pack); }

    /**
     * Convert the given base to its packed value.
     * @param ascii ASCII representation of the base.
     * @return Packed value.
     */
    public static int toPack( byte ascii )
    {
        for( Map.Entry<Integer,Byte> entry: basesByPack.entrySet() ) {
            if( entry.getValue().equals(ascii) )
                return entry.getKey();
        }
        throw new ReviewedStingException(String.format("Base %c is an invalid base to pack", (char)ascii));
    }

    /**
     * Convert the ASCII representation of a base to its 'normalized' representation.
     * @param base The base itself.
     * @return The byte, if present.  Null if unknown.
     */
    public static Byte fromASCII( byte base ) {
        Byte found = null;
        for( Byte normalized: allBases ) {
            if( normalized.equals(base) ) {
                found = normalized;
                break;
            }
        }
        return found;
    }
}
