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
package edu.mit.broad.arachne;


import edu.mit.broad.sam.util.CloseableIterator;

import java.io.*;


/**
 * Reader for arachne Fastb files.
 */
public class FastbReader
    implements CloseableIterator<String> {

    // Notes on fastb file format
    //
    // Fastb files contain the serialized contents of an arachne vecbasevector,
    // which is a typedef for mastervec<basevector, unsigned int>.
    // The serialization of mastervec objects starts with a 24 byte mv_file_control_block,
    // followed by N variable length segments (one for each element of the mastervec vector),
    // followed by an offset table containing N 8-byte file offsets to the N variable length
    // segments, followed by N fixed length data segments, one for each vector element.
    // Thus, reading a single element of the mastervec vector requires reading from three
    // separate places in the file (the offset table, the variable length section and the
    // fixed length section).
    //
    // The mastervec file header is 24 bytes arranged as follows:
    //  n               4-byte signed(?) integer (number of entries)
    //  c1              1-byte unsigned bit mask (see below)
    //  reserved        1-byte unused
    //  sizeX           1-byte unsigned, sizeof first template parameter (16 for fastb files)
    //  sizeA           1-byte unsigned, sizeof second template parameter (4 for fastb files)
    //  offsets_start   8-byte signed(?) integer, file offset of offset table
    //  static_start    8-byte signed(?) integer, file offset of static data (fixed size section)
    //
    // For fastb files, the fixed size section contains 4 bytes for each object, which is the
    // unsigned(?) count of the number of bases in this entry.
    // For fastb files, the variable length section contains a bit vector with two bits per base.
    // The bases are encoded as follows: A = 0, C = 1, G = 2, T = 3.
    //
    // For fastb files, in the file header N is the number of entries in the fastb file.
    // c1 is unused/unimplemented except that the two low-order bits should be 0x01, indicating
    // that we are using the single-file representation.  There is also apparently a three-file
    // representation that looks the same except that the offset table and static (fixed length)
    // table are in separate files named <basename>.offsets and <basename>.static.
    // The sizeX should be 16 for fastb files and sizeA should be 4.
    //
    // Note that in fastb files, the sequences are not identified by name or id, only by index
    // (zero based) into the mastervec object.  There is no representation for bases other than
    // ACGT (i.e. Ns cannot be encoded).

    private static final char[] BASES = { 'A', 'C', 'G', 'T' };

    private File mFile;
    private RandomAccessFile mRandomFile;
    private int mEntryCount;
    private long mOffsetTableOffset;
    private long mLengthTableOffset;
    private int mCurrentPosition;
    private byte[] mIOBuffer = new byte[8];


    public FastbReader(File file)
        throws IOException {
        mFile = file;
        mRandomFile = new RandomAccessFile(mFile, "r");
        readHeader();
    }

    public int getSequenceCount() {
        return mEntryCount;
    }

    public boolean hasNext() {
        return (mCurrentPosition < mEntryCount);
    }

    public String next() {
        if (!hasNext()) {
            throw new IllegalStateException("Iterator exhausted");
        }
        try {
            return readSequence(mCurrentPosition);
        } catch (IOException exc) {
            throw new RuntimeException(exc.getMessage(), exc);
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Not supported: remove");
    }

    public void close() {
        if (mRandomFile != null) {
            mEntryCount = 0;
            mCurrentPosition = 0;
            try {
                mRandomFile.close();
            } catch (IOException exc) {
                throw new RuntimeException(exc.getMessage(), exc);
            } finally {
                mRandomFile = null;
            }
        }
    }

    public String readSequence(int n)
        throws IOException {
        if (mRandomFile == null) {
            throw new IllegalStateException("Reader is closed");
        }
        if (n < 0 || n >= mEntryCount) {
            throw new IndexOutOfBoundsException("Illegal index: " + n);
        }
        long offset = getEntryOffset(n);
        int length = getEntryBaseCount(n);
        String result = readBases(offset, length);
        mCurrentPosition = n+1;
        return result;
    }

    private void readHeader()
        throws IOException {

        byte[] fileControlBlock = new byte[24];
        mRandomFile.readFully(fileControlBlock, 0, 24);

        int word2 = deserializeInt(fileControlBlock, 4);
        int nFiles = word2 & 0x3;
        int sizeX = (word2 >> 16) & 0xFF;
        int sizeA = (word2 >> 24) & 0xFF;
        if (nFiles != 1) {
            throw new RuntimeException(mFile + ": Invalid file header: nFiles = " + nFiles);
        }
        if (sizeX != 16) {
            throw new RuntimeException(mFile + ": Invalid file header: sizeX = " + sizeX);
        }
        if (sizeA != 4) {
            throw new RuntimeException(mFile + ": Invalid file header: sizeX = " + sizeA);
        }
        mEntryCount = deserializeInt(fileControlBlock, 0);
        mOffsetTableOffset = deserializeLong(fileControlBlock, 8);
        mLengthTableOffset = deserializeLong(fileControlBlock, 16);
    }

    private long getEntryOffset(int n)
        throws IOException {
        mRandomFile.seek(mOffsetTableOffset + 8 * n);
        mRandomFile.readFully(mIOBuffer, 0, 8);
        return deserializeLong(mIOBuffer, 0);
    }

    private int getEntryBaseCount(int n)
        throws IOException {
        mRandomFile.seek(mLengthTableOffset + 4 * n);
        mRandomFile.readFully(mIOBuffer, 0, 4);
        return deserializeInt(mIOBuffer, 0);
    }

    private String readBases(long fileOffset, int baseCount)
        throws IOException {


        int byteCount = (baseCount + 3) / 4;
        byte[] data = new byte[byteCount];
        mRandomFile.seek(fileOffset);
        mRandomFile.readFully(data, 0, byteCount);

        int baseIndex = 0;
        int dataIndex = 0;
        char[] baseBuffer = new char[baseCount];
        while (baseIndex < baseCount) {
            int b = data[dataIndex++];
            int count = Math.min(4, baseCount - baseIndex);
            for (int i = 0; i < count; i++) {
                baseBuffer[baseIndex++] = BASES[b & 0x3];
                b = b >> 2;
            }
        }
        return new String(baseBuffer);
    }

    private int deserializeInt(byte[] buffer, int offset) {
        int byte1 = buffer[offset] & 0xFF;
        int byte2 = buffer[offset+1] & 0xFF;
        int byte3 = buffer[offset+2] & 0xFF;
        int byte4 = buffer[offset+3] & 0xFF;
        return (byte1 | (byte2 << 8) | (byte3 << 16) | (byte4 << 24));
    }

    private long deserializeLong(byte[] buffer, int offset) {
        long int1 = deserializeInt(buffer, offset) & 0xFFFFFFFFL;
        long int2 = deserializeInt(buffer, offset+4) & 0xFFFFFFFFL;
        return (int1 | (int2 << 32));
    }

    // Stub for interactive use (see also Fastb2Fasta)
    public static void main(String[] args)
        throws Exception {
        FastbReader reader = new FastbReader(new File(args[0]));
        int readId = 0;
        while (reader.hasNext()) {
            System.out.println(">" + readId);
            System.out.println(reader.next());
            readId++;
        }
        reader.close();
    }
}

