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

package org.broadinstitute.gatk.engine.datasources.reads.utilities;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Calculates the bounds of each BGZF block in a BAM index file, along with
 */
public class PrintBGZFBounds extends CommandLineProgram {
    @Argument(fullName="input",shortName="I",doc="Input bai file to process",required=true)
    private File input = null;

    private final int BYTE_SIZE_IN_BYTES = Byte.SIZE / 8;
    private final int INT_SIZE_IN_BYTES = Integer.SIZE / 8;
    private final int SHORT_SIZE_IN_BYTES = INT_SIZE_IN_BYTES / 2;

    /**
     * ID1 + ID2 + CM + FLG + MTIME + XFL + OS + XLEN.
     */
    private final int HEADER_SIZE = BYTE_SIZE_IN_BYTES*4+INT_SIZE_IN_BYTES+BYTE_SIZE_IN_BYTES*2+SHORT_SIZE_IN_BYTES + BYTE_SIZE_IN_BYTES*2 + SHORT_SIZE_IN_BYTES*2;;

    /**
     * CRC32 + ISIZE
     */
    private final int FOOTER_SIZE = INT_SIZE_IN_BYTES*2;

    @Override
    public int execute() throws IOException {
        FileInputStream fis = new FileInputStream(input);
        ByteBuffer headerBuffer = allocateBuffer(HEADER_SIZE);
        ByteBuffer footerBuffer = allocateBuffer(FOOTER_SIZE);

        float compressedSize = 0;
        float uncompressedSize = 0;
        long totalBlocks = 0;

        //SAMFileReader reader = new SAMFileReader(input);

        while(true) {
            final long blockStart = fis.getChannel().position();

            int totalRead = fis.getChannel().read(headerBuffer);
            if(totalRead <= 0)
                break;
            headerBuffer.flip();

            // Read out header information, including subfield IDs.
            headerBuffer.position(headerBuffer.capacity()-BYTE_SIZE_IN_BYTES*2);
            final int cDataSize = headerBuffer.getShort()-HEADER_SIZE-FOOTER_SIZE+1;
            compressedSize += cDataSize;

            // Skip past body.
            fis.getChannel().position(fis.getChannel().position()+cDataSize);

            // Read the footer
            fis.getChannel().read(footerBuffer);
            footerBuffer.flip();

            // Retrieve the uncompressed size from the footer.
            footerBuffer.position(footerBuffer.capacity()-INT_SIZE_IN_BYTES);
            uncompressedSize += footerBuffer.getInt();

            // Reset buffers for subsequent reads.
            headerBuffer.flip();
            footerBuffer.flip();

            totalBlocks++;

            final long blockStop = fis.getChannel().position() - 1;

            System.out.printf("BGZF block %d: [%d-%d]%n",totalBlocks,blockStart,blockStop);
        }

        System.out.printf("SUCCESS!  Average compressed block size = %f, average uncompressed size = %f, compressed/uncompressed ratio: %f%n",compressedSize/totalBlocks,uncompressedSize/totalBlocks,compressedSize/uncompressedSize);

        return 0;
    }

    private ByteBuffer allocateBuffer(final int size) {
        ByteBuffer buffer = ByteBuffer.allocate(size);
        buffer.order(ByteOrder.LITTLE_ENDIAN);
        return buffer;
    }

    /**
     * Required main method implementation.
     * @param argv Command-line argument text.
     * @throws Exception on error.
     */
    public static void main(String[] argv) throws Exception {
        int returnCode = 0;
        try {
            PrintBGZFBounds instance = new PrintBGZFBounds();
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
