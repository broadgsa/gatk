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

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * A general-purpose stream for writing packed bases.
 *
 * @author mhanna
 * @version 0.1
 */
public class BasePackedOutputStream<T> {
    /**
     * Type of object to pack.
     */
    private final Class<T> type;

    /**
     * How many bases can be stored in the given data structure?
     */
    private final int basesPerType;

    /**
     * Ultimate target for the packed bases.
     */
    private final OutputStream targetOutputStream;

    /**
     * A fixed-size buffer for word-packed data.
     */
    private final ByteBuffer buffer;

    public BasePackedOutputStream( Class<T> type, File outputFile, ByteOrder byteOrder ) throws FileNotFoundException {
        this(type,new BufferedOutputStream(new FileOutputStream(outputFile)),byteOrder);
    }

    /**
     * Write packed bases to the given output stream.
     * @param type Type of data to pack bases into.
     * @param outputStream Output stream to which to write packed bases.
     * @param byteOrder Switch between big endian / little endian when reading / writing files.
     */
    public BasePackedOutputStream( Class<T> type, OutputStream outputStream, ByteOrder byteOrder) {
        this.targetOutputStream = outputStream;
        this.type = type;
        basesPerType = PackUtils.bitsInType(type)/PackUtils.BITS_PER_BASE;
        this.buffer = ByteBuffer.allocate(basesPerType/PackUtils.ALPHABET_SIZE).order(byteOrder);
    }

    /**
     * Writes the given base to the output stream.  Will write only this base; no packing will be performed.
     * @param base List of bases to write.
     * @throws IOException if an I/O error occurs.
     */
    public void write( int base ) throws IOException {
        write( new byte[] { (byte)base } );
    }

    /**
     * Writes an array of bases to the target output stream.
     * @param bases List of bases to write.
     * @throws IOException if an I/O error occurs.
     */
    public void write( byte[] bases ) throws IOException {
        write(bases,0,bases.length);
    }

    /**
     * Writes a subset of the array of bases to the output stream.
     * @param bases List of bases to write.
     * @param offset site at which to start writing.
     * @param length number of bases to write.
     * @throws IOException if an I/O error occurs.
     */
    public void write( byte[] bases, int offset, int length ) throws IOException {
        int packedBases = 0;
        int positionInPack = 0;

        for( int base = offset; base < offset+length; base++ ) {
            packedBases = packBase(bases[base], packedBases, positionInPack);

            // Increment the packed counter.  If all possible bases have been squeezed into this byte, write it out.
            positionInPack = ++positionInPack % basesPerType;
            if( positionInPack == 0 ) {
                writePackedBases(packedBases);
                packedBases = 0;
            }
        }

        if( positionInPack > 0 )
            writePackedBases(packedBases);
    }

    /**
     * Flush the contents of the OutputStream to disk.
     * @throws IOException if an I/O error occurs.
     */
    public void flush() throws IOException {
        targetOutputStream.flush();
    }

    /**
     * Closes the given output stream.
     * @throws IOException if an I/O error occurs.
     */
    public void close() throws IOException {
        targetOutputStream.close();
    }

    /**
     * Pack the given base into the basepack.
     * @param base The base to pack.
     * @param basePack Target for the pack operation.
     * @param position Position within the pack to which to add the base.
     * @return The packed integer.
     */
    private int packBase( byte base, int basePack, int position ) {
        basePack |= (PackUtils.packBase(base) << 2*(basesPerType-position-1));
        return basePack;
    }    

    /**
     * Write the given packed base structure to the output file.
     * @param packedBases Packed bases to write.
     * @throws IOException on error writing to the file.
     */
    private void writePackedBases(int packedBases) throws IOException {
        buffer.rewind();
        if( type == Integer.class )
            buffer.putInt(packedBases);
        else if( type == Byte.class )
            buffer.put((byte)packedBases);
        else
            throw new ReviewedGATKException("Cannot pack bases into type " + type.getName());
        targetOutputStream.write(buffer.array());        
    }
}
