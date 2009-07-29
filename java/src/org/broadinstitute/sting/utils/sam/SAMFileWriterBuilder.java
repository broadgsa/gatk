/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

import java.io.File;

import org.broadinstitute.sting.utils.StingException;

/**
 * Allows the user to steadily accumulate information about what
 * components go into a SAM file writer, ultimately using this
 * information to create a SAM file writer on demand.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMFileWriterBuilder {
    /**
     * Default compression level for newly constructed SAM files.
     * Default to 5 (based on research by Alec Wysoker)
     */
    public static final int DEFAULT_COMPRESSION_LEVEL = 5;

    /**
     * To which file should output be written?
     */
    private File samFile = null;

    /**
     * Which header should be used when writing the SAM file?
     */
    private SAMFileHeader header = null;

    /**
     * What compression level should be used when building this file?
     */
    private int compressionLevel = DEFAULT_COMPRESSION_LEVEL;

    /**
     * Sets the handle of the sam file to which data should be written.
     * @param samFile The SAM file into which data should flow.
     */
    public void setSAMFile( File samFile ) {
        this.samFile = samFile;
    }

    /**
     * Sets the header to be written at the head of this SAM file.
     * @param header Header to write.
     */
    public void setSAMFileHeader( SAMFileHeader header ) {
        this.header = header;
    }

    /**
     * Sets the compression level to use when writing this BAM file.
     * @param compressionLevel Compression level to use when writing this SAM file.
     */
    public void setCompressionLevel( int compressionLevel ) {
        this.compressionLevel = compressionLevel;
    }

    /**
     * Create the SAM writer, given the constituent parts accrued.
     * @return Newly minted SAM file writer.
     */
    public SAMFileWriter build() {
        if( samFile == null )
            throw new StingException( "Filename for output sam file must be supplied.");
        if( header == null )
            throw new StingException( "Header for output sam file must be supplied.");
        return new SAMFileWriterFactory().makeBAMWriter( header, true, samFile, compressionLevel );
    }
}
