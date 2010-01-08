package net.sf.samtools;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

/**
 * Test harness for playing with sharding by BGZF block.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockTestHarness {
    private static void usage() {
        System.out.println("Usage: ChunkTestHarness <filename>.bam");
        System.exit(1);
    }

    public static void main(String args[]) throws IOException {
        if(args.length == 0)
            usage();

        String bamFileName = args[0];
        if(!bamFileName.endsWith(".bam"))
            usage();

        File bamFile = new File(bamFileName);
        if(!bamFile.exists())
            usage();

        Chunk headerLocation = BAMFileHeaderLoader.load(bamFile).getLocation();
        System.out.printf("Header location = %s%n", headerLocation);

        BAMChunkIterator chunkIterator = new BAMChunkIterator(new BAMBlockIterator(bamFile),
                                                              Collections.singletonList(headerLocation));
        while(chunkIterator.hasNext()) {
            Chunk chunk = chunkIterator.next();
        }

        // The following test code blocks a bunch of files together.
        /*
        BAMBlockIterator blockIterator = new BAMBlockIterator(bamFile);
        long blockCount = 0;

        long startTime = System.currentTimeMillis();
        while(blockIterator.hasNext()) {
            Block block = blockIterator.next();
            blockCount++;
            System.out.println(block);
        }
        long endTime = System.currentTimeMillis();

        System.out.printf("Number of chunks: %d; elapsed time: %dms%n", blockCount, endTime-startTime);
        */
    }
}
