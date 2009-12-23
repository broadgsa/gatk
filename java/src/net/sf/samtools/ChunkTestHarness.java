package net.sf.samtools;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;

/**
 * Test harness for playing with sharding by BGZF block.
 *
 * @author mhanna
 * @version 0.1
 */
public class ChunkTestHarness {
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

        FileInputStream bamInputStream = new FileInputStream(bamFile);
        FileChannel bamInputChannel = bamInputStream.getChannel();

        BAMChunkIterator chunkIterator = new BAMChunkIterator(bamInputChannel);
        long chunkCount = 0;

        long startTime = System.currentTimeMillis();
        while(chunkIterator.hasNext()) {
            Chunk chunk = chunkIterator.next();
            chunkCount++;
            System.out.printf("Chunk: [%d,%d)\tByte offsets: [%d,%d)%n",chunk.getChunkStart(),
                                                                        chunk.getChunkEnd(),
                                                                        chunk.getChunkStart()>>16,
                                                                        chunk.getChunkEnd()>>16);
        }
        long endTime = System.currentTimeMillis();

        System.out.printf("Number of chunks: %d; elasped time: %dms%n", chunkCount, endTime-startTime);
    }
}
