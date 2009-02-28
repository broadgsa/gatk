/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.util;

import edu.mit.broad.sam.util.SortingCollection;
import edu.mit.broad.sam.util.RuntimeIOException;

import java.util.Comparator;
import java.nio.ByteBuffer;
import java.io.OutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.File;

/**
 * Factory to create new String SortingCollections
 *
 * @author Kathleen Tibbetts
 */
public class StringSortingCollectionFactory {

    private static final File TEMP_DIR = new File(System.getProperty("java.io.tmpdir"), "StringSortingCollectionFactory");
    private static final int MAX_RECORDS_IN_RAM = 20000;
    
    private StringSortingCollectionFactory() {
    }

    public static SortingCollection<String> newCollection() {
        return SortingCollection.newInstance(
                String.class, new StringCodec(), new StringComparator(), MAX_RECORDS_IN_RAM, TEMP_DIR);
    }

    static class StringCodec implements SortingCollection.Codec<String> {
        ByteBuffer byteBuffer = ByteBuffer.allocate(4);
        OutputStream os;
        InputStream is;

        /** Returns a new StringCodec. */
        public SortingCollection.Codec<String> clone() {
            return new StringCodec();
        }

        /**
         * Where to write encoded output
         *
         * @param os    the output stream to encode output
         */
        public void setOutputStream(final OutputStream os) {
            this.os = os;
        }

        /**
         * Where to read encoded input from
         *
         * @param is where to read encoded input from
         */
        public void setInputStream(final InputStream is) {
            this.is = is;
        }

        /**
         * Write object to file
         *
         * @param val what to write
         */
        public void encode(final String val) {
            try {
                byteBuffer.clear();
                byteBuffer.putInt(val.length());
                os.write(byteBuffer.array());
                os.write(val.getBytes());
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        /**
         * Read the next record from the input stream and convert into a java object.
         *
         * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
         *         a record.
         */
        public String decode() {
            try {
                byteBuffer.clear();
                int bytesRead = is.read(byteBuffer.array());
                if (bytesRead == -1) {
                    return null;
                }
                if (bytesRead != 4) {
                    throw new RuntimeException("Unexpected EOF in middle of record");
                }
                byteBuffer.limit(4);
                final int length = byteBuffer.getInt();
                final byte[] buf = new byte[length];
                bytesRead = is.read(buf);
                if (bytesRead != length) {
                    throw new RuntimeException("Unexpected EOF in middle of record");
                }
                return new String(buf);
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }
    }

    static class StringComparator implements Comparator<String> {

        public int compare(final String s, final String s1) {
            return s.compareTo(s1);
        }
    }

}
