/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.bcf2;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broad.tribble.FeatureCodec;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

public class BCF2Decoder {
    final protected static Logger logger = Logger.getLogger(FeatureCodec.class);

    byte[] recordBytes;
    ByteArrayInputStream recordStream;

    public BCF2Decoder() {
        // nothing to do
    }

    /**
     * Create a new decoder ready to read BCF2 data from the byte[] recordBytes, for testing purposes
     *
     * @param recordBytes
     */
    protected BCF2Decoder(final byte[] recordBytes) {
        setRecordBytes(recordBytes);
    }

    // ----------------------------------------------------------------------
    //
    // Routines to load, set, skip blocks of underlying data we are decoding
    //
    // ----------------------------------------------------------------------

    /**
     * Reads the next record from input stream and prepare this decoder to decode values from it
     *
     * @param stream
     * @return
     */
    public void readNextBlock(final int blockSizeInBytes, final InputStream stream) {
        setRecordBytes(readRecordBytes(blockSizeInBytes, stream));
    }

    /**
     * Skips the next record from input stream, invalidating current block data
     *
     * @param stream
     * @return
     */
    public void skipNextBlock(final int blockSizeInBytes, final InputStream stream) {
        try {
            final int bytesRead = (int)stream.skip(blockSizeInBytes);
            validateReadBytes(bytesRead, blockSizeInBytes);
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 file", e);
        }
        this.recordBytes = null;
        this.recordStream = null;
    }

    /**
     * Returns the byte[] for the block of data we are currently decoding
     * @return
     */
    public byte[] getRecordBytes() {
        return recordBytes;
    }

    /**
     * The size of the current block in bytes
     *
     * @return
     */
    public int getBlockSize() {
        return recordBytes.length;
    }

    public boolean blockIsFullyDecoded() {
        return recordStream.available() == 0;
    }

    /**
     * Use the recordBytes[] to read BCF2 records from now on
     *
     * @param recordBytes
     */
    public void setRecordBytes(final byte[] recordBytes) {
        assert recordBytes != null;

        this.recordBytes = recordBytes;
        this.recordStream = new ByteArrayInputStream(recordBytes);
    }

    // ----------------------------------------------------------------------
    //
    // High-level decoder
    //
    // ----------------------------------------------------------------------

    public final Object decodeTypedValue() {
        final byte typeDescriptor = readTypeDescriptor();
        return decodeTypedValue(typeDescriptor);
    }

    public final Object decodeTypedValue(final byte typeDescriptor) {
        final int size = decodeNumberOfElements(typeDescriptor);
        final BCF2Type type = BCF2Utils.decodeType(typeDescriptor);

        assert size >= 0;

        if ( size == 0 ) {
            return null;
        } else if ( type == BCF2Type.CHAR ) { // special case string decoding for efficiency
            return decodeLiteralString(size);
        } else if ( size == 1 ) {
            return decodeSingleValue(type);
        } else {
            final ArrayList<Object> ints = new ArrayList<Object>(size);
            for ( int i = 0; i < size; i++ ) {
                final Object val = decodeSingleValue(type);
                if ( val == null ) continue; // auto-pruning.  We remove trailing nulls
                ints.add(val);
            }
            return ints.isEmpty() ? null : ints; // return null when all of the values are null
        }
    }

    public final Object decodeSingleValue(final BCF2Type type) {
        // TODO -- decodeTypedValue should integrate this routine
        final int value = decodeInt(type);

        if ( value == type.getMissingBytes() )
            return null;
        else {
            switch (type) {
                case INT8:
                case INT16:
                case INT32: return value;
                case FLOAT: return rawFloatToFloat(value);
                case CHAR:  return value & 0xFF; // TODO -- I cannot imagine why we'd get here, as string needs to be special cased
                default:    throw new ReviewedStingException("BCF2 codec doesn't know how to decode type " + type );
            }
        }
    }

    // ----------------------------------------------------------------------
    //
    // Decode raw primitive data types (ints, floats, and strings)
    //
    // ----------------------------------------------------------------------

    private final Object decodeLiteralString(final int size) {
        assert size > 0;

        // TODO -- assumes size > 0
        final byte[] bytes = new byte[size]; // TODO -- in principle should just grab bytes from underlying array
        try {
            recordStream.read(bytes);
            final String s = new String(bytes);
            return BCF2Utils.isCollapsedString(s) ? BCF2Utils.exploreStringList(s) : s;
        } catch ( IOException e ) {
            throw new ReviewedStingException("readByte failure", e);
        }
    }

    @Ensures("result >= 0")
    public final int decodeNumberOfElements(final byte typeDescriptor) {
        if ( BCF2Utils.sizeIsOverflow(typeDescriptor) )
            // -1 ensures we explode immediately with a bad size if the result is missing
            return decodeInt(readTypeDescriptor(), -1);
        else
            // the size is inline, so just decode it
            return BCF2Utils.decodeSize(typeDescriptor);
    }

    /**
     * Decode an int from the stream.  If the value in the stream is missing,
     * returns missingValue.  Requires the typeDescriptor indicate an inline
     * single element event
     *
     * @param typeDescriptor
     * @return
     */
    @Requires("BCF2Utils.decodeSize(typeDescriptor) == 1")
    public final int decodeInt(final byte typeDescriptor, final int missingValue) {
        final BCF2Type type = BCF2Utils.decodeType(typeDescriptor);
        final int i = decodeInt(type);
        return i == type.getMissingBytes() ? missingValue : i;
    }

    @Requires("type != null")
    public final int decodeInt(final BCF2Type type) {
        return BCF2Utils.readInt(type.getSizeInBytes(), recordStream);
    }

    /**
     * Low-level reader for int[]
     *
     * Requires a typeDescriptor so the function knows how many elements to read,
     * and how they are encoded.
     *
     * If size == 0 => result is null
     * If size > 0 => result depends on the actual values in the stream
     *      -- If the first element read is MISSING, result is null (all values are missing)
     *      -- Else result = int[N] where N is the first N non-missing values decoded
     *
     * @param maybeDest if not null we'll not allocate space for the vector, but instead use
     *                  the externally allocated array of ints to store values.  If the
     *                  size of this vector is < the actual size of the elements, we'll be
     *                  forced to use freshly allocated arrays.  Also note that padded
     *                  int elements are still forced to do a fresh allocation as well.
     * @return see description
     */
    @Requires({"BCF2Type.INTEGERS.contains(type)", "size >= 0"})
    public final int[] decodeIntArray(final int size, final BCF2Type type, int[] maybeDest) {
        if ( size == 0 ) {
            return null;
        } else {
            if ( maybeDest != null && maybeDest.length < size )
                maybeDest = null; // by nulling this out we ensure that we do fresh allocations as maybeDest is too small

            final int val1 = decodeInt(type);
            if ( val1 == type.getMissingBytes() ) {
                // fast path for first element being missing
                for ( int i = 1; i < size; i++ ) decodeInt(type);
                return null;
            } else {
                // we know we will have at least 1 element, so making the int[] is worth it
                final int[] ints = maybeDest == null ? new int[size] : maybeDest;
                ints[0] = val1; // we already read the first one
                for ( int i = 1; i < size; i++ ) {
                    ints[i] = decodeInt(type);
                    if ( ints[i] == type.getMissingBytes() ) {
                        // read the rest of the missing values, dropping them
                        for ( int j = i + 1; j < size; j++ ) decodeInt(type);
                        // deal with auto-pruning by returning an int[] containing
                        // only the non-MISSING values.  We do this by copying the first
                        // i elements, as i itself is missing
                        return Arrays.copyOf(ints, i);
                    }
                }
                return ints; // all of the elements were non-MISSING
            }
        }
    }

    public final int[] decodeIntArray(final byte typeDescriptor) {
        final int size = decodeNumberOfElements(typeDescriptor);
        final BCF2Type type = BCF2Utils.decodeType(typeDescriptor);
        return decodeIntArray(size, type, null);
    }

        public final double rawFloatToFloat(final int rawFloat) {
        return (double)Float.intBitsToFloat(rawFloat);
    }

    // ----------------------------------------------------------------------
    //
    // Utility functions
    //
    // ----------------------------------------------------------------------

    /**
     * Read the size of the next block from inputStream
     *
     * @param inputStream
     * @return
     */
    public final int readBlockSize(final InputStream inputStream) {
        return BCF2Utils.readInt(4, inputStream);
    }

    /**
     *
     * @param inputStream
     * @return
     */
    private final static byte[] readRecordBytes(final int blockSizeInBytes, final InputStream inputStream) {
        assert blockSizeInBytes >= 0;

        final byte[] record = new byte[blockSizeInBytes];
        try {
            final int bytesRead = inputStream.read(record);
            validateReadBytes(bytesRead, blockSizeInBytes);
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 file", e);
        }

        return record;
    }

    private final static void validateReadBytes(final int actuallyRead, final int expected) {
        assert expected >= 0;

        if ( actuallyRead < expected ) {
            throw new UserException.MalformedBCF2(String.format("Failed to read next complete record: %s",
                    actuallyRead == -1 ?
                            "premature end of input stream" :
                            String.format("expected %d bytes but read only %d", expected, actuallyRead)));
        }
    }

    public final byte readTypeDescriptor() {
        return BCF2Utils.readByte(recordStream);
    }
}
