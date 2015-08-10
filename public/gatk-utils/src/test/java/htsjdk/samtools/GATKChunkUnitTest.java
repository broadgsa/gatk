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

package htsjdk.samtools;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Test basic functionality of the GATK chunk, giving informative size capabilities, etc.
 */
public class GATKChunkUnitTest {
    private static final int FULL_BLOCK_COMPRESSED_SIZE = 25559;
    private static final int FULL_BLOCK_UNCOMPRESSED_SIZE = 65536;
    private static final int HALF_BLOCK_UNCOMPRESSED_SIZE = FULL_BLOCK_UNCOMPRESSED_SIZE/2;

    @Test
    public void testSizeOfEmptyChunk() {
        GATKChunk chunk = new GATKChunk(0,0);
        Assert.assertEquals(chunk.size(),0,"Empty chunk's size is not equal to 0.");
    }

    @Test
    public void testSizeOfChunkWithinSingleBlock() {
        GATKChunk chunk = new GATKChunk(0,FULL_BLOCK_UNCOMPRESSED_SIZE-1);
        Assert.assertEquals(chunk.size(),FULL_BLOCK_UNCOMPRESSED_SIZE-1,"Chunk spanning limits of block is returning wrong size.");

        chunk = new GATKChunk(0,HALF_BLOCK_UNCOMPRESSED_SIZE);
        Assert.assertEquals(chunk.size(),HALF_BLOCK_UNCOMPRESSED_SIZE,"Chunk spanning 1/2 block is returning the wrong size.");
    }

    @Test
    public void testSizeOfSingleBlock() {
        GATKChunk chunk = new GATKChunk(0,FULL_BLOCK_COMPRESSED_SIZE<<16);
        Assert.assertEquals(chunk.size(),FULL_BLOCK_UNCOMPRESSED_SIZE,"Chunk spanning complete block returns incorrect size.");
    }

    @Test
    public void testSizeOfBlockAndAHalf() {
        GATKChunk chunk = new GATKChunk(0,(FULL_BLOCK_COMPRESSED_SIZE<<16)+HALF_BLOCK_UNCOMPRESSED_SIZE);
        Assert.assertEquals(chunk.size(),FULL_BLOCK_UNCOMPRESSED_SIZE+HALF_BLOCK_UNCOMPRESSED_SIZE,"Chunk spanning 1.5 blocks returns incorrect size.");
    }

    @Test
    public void testSizeOfHalfBlock() {
        GATKChunk chunk = new GATKChunk(HALF_BLOCK_UNCOMPRESSED_SIZE,FULL_BLOCK_COMPRESSED_SIZE<<16);
        Assert.assertEquals(chunk.size(),HALF_BLOCK_UNCOMPRESSED_SIZE,"Chunk spanning 0.5 blocks returns incorrect size.");
    }
}
