/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.alignment.reference.bwt;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.*;

/**
 * Enhanced enum representation of a base.
 *
 * @author mhanna
 * @version 0.1
 */
public class Bases implements Iterable<Byte>
{
    public static final byte A = 'A';
    public static final byte C = 'C';
    public static final byte G = 'G';
    public static final byte T = 'T';

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
        throw new ReviewedGATKException(String.format("Base %c is an invalid base to pack", (char)ascii));
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
