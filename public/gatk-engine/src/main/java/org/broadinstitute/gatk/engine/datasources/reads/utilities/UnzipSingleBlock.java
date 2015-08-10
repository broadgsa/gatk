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

import htsjdk.samtools.util.BlockGunzipper;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.commandline.Input;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Test decompression of a single BGZF block.
 */
public class UnzipSingleBlock extends CommandLineProgram {
    @Input(fullName = "block_file", shortName = "b", doc = "block file over which to test unzipping", required = true)
    private File blockFile;

    @Input(fullName = "compressed_block_size", shortName = "cbs", doc = "size of compressed block", required = true)
    private int compressedBufferSize;

    public int execute() throws IOException, NoSuchMethodException, IllegalAccessException, InvocationTargetException {
        byte[] compressedBuffer = new byte[(int)blockFile.length()];
        byte[] uncompressedBuffer = new byte[65536];

        FileInputStream fis = new FileInputStream(blockFile);
        fis.read(compressedBuffer);
        fis.close();

        BlockGunzipper gunzipper = new BlockGunzipper();
        gunzipper.setCheckCrcs(true);
        Method unzipBlock = BlockGunzipper.class.getDeclaredMethod("unzipBlock",byte[].class,byte[].class,Integer.TYPE);
        unzipBlock.setAccessible(true);

        unzipBlock.invoke(gunzipper,uncompressedBuffer,compressedBuffer,compressedBufferSize);

        System.out.printf("SUCCESS!%n");

        return 0;
    }

    /**
     * Required main method implementation.
     * @param argv Command-line argument text.
     * @throws Exception on error.
     */
    public static void main(String[] argv) throws Exception {
        int returnCode = 0;
        try {
            UnzipSingleBlock instance = new UnzipSingleBlock();
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
