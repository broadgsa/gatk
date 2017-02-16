/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.datasources.reads.utilities;

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.commandline.Input;

import java.io.File;
import java.io.IOException;

/**
 * Test decompression of a single BGZF block.
 */
public class UnzipSingleBlock extends CommandLineProgram {
    @Input(fullName = "block_file", shortName = "b", doc = "block file over which to test unzipping", required = true)
    private File blockFile;

    @Input(fullName = "compressed_block_size", shortName = "cbs", doc = "size of compressed block", required = true)
    private int compressedBufferSize;

    public int execute() throws IOException {
        final byte[] uncompressedBuffer = new byte[65536];

        final BlockCompressedInputStream gunzipper = new BlockCompressedInputStream(blockFile);
        gunzipper.setCheckCrcs(true);
        gunzipper.read(uncompressedBuffer);

        System.out.printf("SUCCESS!%n");

        return 0;
    }

    /**
     * Required main method implementation.
     * @param argv Command-line argument text.
     * @throws Exception on error.
     */
    public static void main(final String[] argv) throws Exception {
        int returnCode = 0;
        try {
            final UnzipSingleBlock instance = new UnzipSingleBlock();
            start(instance, argv);
            returnCode = 0;
        }
        catch(Exception ex) {
            returnCode = 1;
            ex.printStackTrace();
            throw ex;
        }
        finally {
            System.exit(returnCode);
        }
    }
}
