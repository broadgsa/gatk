package org.broadinstitute.sting.gatk.datasources;

import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.PrintStream;
import java.nio.channels.FileChannel;

import net.sf.samtools.*;

/**
 *
 *
 * @author mhanna
 * @version 0.1
 */
public class BAMFileStat extends CommandLineProgram {
    public enum CommandType { ShowBlocks, ShowIndex }

    @Argument(doc="Which operation to run.",required=true)
    private CommandType command;

    @Argument(doc="The BAM file to inspect.",required=true)
    private File bamFile;

    @Argument(doc="The range of blocks to inspect.",required=false)
    private String range;

    public int execute() {
        Integer startPosition = null, stopPosition = null;
        if(range != null) {
            int dashPosition = range.indexOf('-');
            if(dashPosition > 0) {
                if(dashPosition > 0)
                    startPosition = Integer.valueOf(range.substring(0,dashPosition));
                if(dashPosition < range.length()-1)
                    stopPosition = Integer.valueOf(range.substring(dashPosition+1));
            }
            else
                startPosition = Integer.valueOf(range);
        }

        switch(command) {
            case ShowBlocks:
                showBlocks(bamFile,startPosition,stopPosition);
                break;
            case ShowIndex:
                showIndexBins(bamFile);
                break;
        }
        return 0;
    }

    /**
     * Required main method implementation.
     * @param argv Command-line arguments.
     */
    public static void main(String[] argv) {
        try {
            BAMFileStat instance = new BAMFileStat();
            start(instance, argv);
            System.exit(CommandLineProgram.result);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }

    private void showBlocks(File bamFile, Integer startPosition, Integer stopPosition) {
        int blockNumber = 0;

        try {
            BAMBlockIterator iterator = new BAMBlockIterator(bamFile);
            while(iterator.hasNext()) {
                Block block = iterator.next();
                blockNumber++;

                if(startPosition != null && startPosition > blockNumber) continue;
                if(stopPosition != null && stopPosition < blockNumber) break;

                System.out.printf("Block number = %d; position = %d; compressed size = %d; uncompressed size = %d%n", blockNumber, block.position, block.compressedBlockSize, block.uncompressedBlockSize);
            }
        }
        catch(IOException ex) {
            throw new StingException("Unable to open BAM file");
        }
    }

    private void showIndexBins(File bamFile) {
        BAMFileIndexContentInspector inspector = new BAMFileIndexContentInspector(bamFile);
        inspector.inspect(System.out,null,null);
    }

    private class BAMFileIndexContentInspector extends BAMFileIndex2 {
        public BAMFileIndexContentInspector(File bamFileIndex) {
            super(bamFileIndex);    
        }

        public void inspect(PrintStream outputStream, Integer startPosition, Integer stopPosition) {
            outputStream.printf("Number of reference sequences: %d%n", this.referenceToBins.size());
            for(int referenceSequence: referenceToBins.keySet()) {
                Bin[] bins = referenceToBins.get(referenceSequence);
                outputStream.printf("Reference sequence: %d%n",referenceSequence);
                outputStream.printf("number of bins: %d%n",bins.length);
                for(Bin bin: bins) {
                    outputStream.printf("\tBin: %d, number of chunks: %d%n", bin.binNumber, bin.chunks.size());
                    for(Chunk chunk: bin.chunks)
                        outputStream.printf("\t\tChunk: %s%n", chunk);
                }
                LinearIndex linearIndex = referenceToLinearIndices.get(referenceSequence);
                outputStream.printf("\t\tIndex entries: %d", linearIndex.indexEntries.length);
                for(long indexEntry: linearIndex.indexEntries)
                    outputStream.printf("%d,",indexEntry);
                outputStream.printf("%n");
            }
        }
    }
}
