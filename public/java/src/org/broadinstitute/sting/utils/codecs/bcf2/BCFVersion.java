package org.broadinstitute.sting.utils.codecs.bcf2;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;

/**
 * Simple holder for BCF version information
 *
 * User: depristo
 * Date: 8/2/12
 * Time: 2:16 PM
 */
public class BCFVersion {
    /**
     * BCF2 begins with the MAGIC info BCF_M_m where M is the major version (currently 2)
     * and m is the minor version, currently 1
     */
    public static final byte[] MAGIC_HEADER_START = "BCF".getBytes();

    final int majorVersion;
    final int minorVersion;

    public BCFVersion(int majorVersion, int minorVersion) {
        this.majorVersion = majorVersion;
        this.minorVersion = minorVersion;
    }

    /**
     * @return the major version number of this BCF file
     */
    public int getMajorVersion() {
        return majorVersion;
    }

    /**
     * @return the minor version number of this BCF file
     */
    public int getMinorVersion() {
        return minorVersion;
    }

    /**
     * Return a new BCFVersion object describing the major and minor version of the BCF file in stream
     *
     * Note that stream must be at the very start of the file.
     *
     * @param stream
     * @return a BCFVersion object, or null if stream doesn't contain a BCF file
     * @throws IOException
     */
    public static BCFVersion readBCFVersion(final InputStream stream) throws IOException {
        final byte[] magicBytes = new byte[MAGIC_HEADER_START.length];
        stream.read(magicBytes);
        if ( Arrays.equals(magicBytes, MAGIC_HEADER_START) ) {
            // we're a BCF file
            final int majorByte = stream.read();
            final int minorByte = stream.read();
            return new BCFVersion( majorByte, minorByte );
        } else
            return null;
    }

    /**
     * Write out the BCF magic information indicating this is a BCF file with corresponding major and minor versions
     * @param out
     * @throws IOException
     */
    public void write(final OutputStream out) throws IOException {
        out.write(MAGIC_HEADER_START);
        out.write(getMajorVersion() & 0xFF);
        out.write(getMinorVersion() & 0xFF);
    }

    @Override
    public String toString() {
        return String.format("BCF%d.%d", getMajorVersion(), getMinorVersion());
    }
}
