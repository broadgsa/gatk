/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.utils.gcf;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.codecs.vcf.IndexingVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.*;

/**
 * GCFWriter implementing the VCFWriter interface
 * @author Your Name
 * @since Date created
 */
public class GCFWriter extends IndexingVCFWriter {
    final boolean skipGenotypes;
    final FileOutputStream fileOutputStream;
    final DataOutputStream dataOutputStream;
    final GCFHeaderBuilder gcfHeaderBuilder;
    int nbytes = 0;
    VCFHeader header = null;
    File location;

    // --------------------------------------------------------------------------------
    //
    // Constructors
    //
    // --------------------------------------------------------------------------------

    public GCFWriter(final File location, final SAMSequenceDictionary refDict, boolean enableOnTheFlyIndexing, boolean doNotWriteGenotypes) {
        super(writerName(location, null), location, null, refDict, enableOnTheFlyIndexing);
        this.location = location;
        this.skipGenotypes = doNotWriteGenotypes;

        // write the output
        try {
            fileOutputStream = new FileOutputStream(location);
            dataOutputStream = createDataOutputStream(fileOutputStream);
            gcfHeaderBuilder = new GCFHeaderBuilder();
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(location, e);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // VCFWriter interface functions
    //
    // --------------------------------------------------------------------------------

    @Override
    public void writeHeader(VCFHeader header) {
        this.header = header;
        try {
            nbytes += GCFHeader.writeHeader(dataOutputStream);
        } catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(getStreamName(), "Couldn't write header", e);
        }
    }

    @Override
    public void add(VariantContext vc) {
        super.add(vc);
        GCF gcf = new GCF(gcfHeaderBuilder, vc, skipGenotypes);
        try {
            nbytes += gcf.write(dataOutputStream);
        } catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile(getStreamName(), "Failed to add gcf record " + gcf + " to stream " + getStreamName(), e);
        }
    }

    @Override
    public void close() {
        // todo -- write out VCF header lines
        GCFHeader gcfHeader = gcfHeaderBuilder.createHeader();
        try {
            long headerPosition = nbytes;
            nbytes += gcfHeader.writeFooter(dataOutputStream);
            dataOutputStream.close();
            //System.out.println("Writing forward reference to " + headerPosition);

            RandomAccessFile raFile = new RandomAccessFile(location, "rw");
            raFile.seek(GCFHeader.HEADER_FORWARD_REFERENCE_OFFSET);
            raFile.writeLong(headerPosition);
            raFile.close();
        } catch ( IOException e ) {
            throw new ReviewedStingException("Failed to close GCFWriter " + getStreamName(), e);
        }

        super.close();
    }

    private static final DataOutputStream createDataOutputStream(final OutputStream stream) {
        return new DataOutputStream(new BufferedOutputStream(stream, GCF.BUFFER_SIZE));
    }

}
