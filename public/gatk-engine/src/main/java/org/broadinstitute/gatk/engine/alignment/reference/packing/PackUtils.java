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

package org.broadinstitute.gatk.engine.alignment.reference.packing;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteOrder;

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
     * Writes a reference sequence to a PAC file.
     * @param outputFile Filename for the PAC file.
     * @param referenceSequence Reference sequence to write.
     * @throws IOException If there's a problem writing to the output file.
     */
    public static void writeReferenceSequence( File outputFile, byte[] referenceSequence ) throws IOException {
        OutputStream outputStream = new FileOutputStream(outputFile);

        BasePackedOutputStream<Byte> basePackedOutputStream = new BasePackedOutputStream<Byte>(Byte.class, outputStream, ByteOrder.BIG_ENDIAN);
        basePackedOutputStream.write(referenceSequence);

        outputStream.write(referenceSequence.length%PackUtils.ALPHABET_SIZE);

        outputStream.close();
    }


    /**
     * How many bits can a given type hold?
     * @param type Type to test.
     * @return Number of bits that the given type can hold.
     */
    public static int bitsInType( Class<?> type ) {
        try {
            long typeSize = type.getField("MAX_VALUE").getLong(null) - type.getField("MIN_VALUE").getLong(null)+1;
            long intTypeSize = (long)Integer.MAX_VALUE - (long)Integer.MIN_VALUE + 1;
            if( typeSize > intTypeSize )
                throw new ReviewedGATKException("Cannot determine number of bits available in type: " + type.getName());
            return (int)(Math.log(typeSize)/Math.log(2));
        }
        catch( NoSuchFieldException ex ) {
            throw new ReviewedGATKException("Cannot determine number of bits available in type: " + type.getName(),ex);
        }
        catch( IllegalAccessException ex ) {
            throw new ReviewedGATKException("Cannot determine number of bits available in type: " + type.getName(),ex);
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
                throw new ReviewedGATKException("Unknown base type: " + base);
        }
    }

    /**
     * Converts a two-bit representation of a base into an ASCII representation of a base. 
     * @param pack Byte from 0-3 indicating which base is represented.
     * @return An ASCII value representing the packed base.
     */
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
                throw new ReviewedGATKException("Unknown pack type: " + pack);
        }
    }

    /**
     * Reverses an unpacked sequence of bases.
     * @param bases bases to reverse.
     */
    public static void reverse( byte[] bases ) {
        for( int i = 0, j = bases.length-1; i < j; i++, j-- ) {
            byte temp = bases[j];
            bases[j] = bases[i];
            bases[i] = temp;
        }        
    }

    /**
     * Given a structure of size <code>size</code> that should be split
     * into <code>partitionSize</code> partitions, how many partitions should
     * be created?  Size of last partition will be <= partitionSize.
     * @param size Total size of the data structure.
     * @param partitionSize Size of an individual partition.
     * @return Number of partitions that would be created.
     */
    public static int numberOfPartitions( long size, long partitionSize ) {
        return (int)((size+partitionSize-1) / partitionSize);    
    }
}
