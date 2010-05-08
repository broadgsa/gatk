/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.datasources;

import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.PrintStream;

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
    private String bamFileName;

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
                throw new StingException("The BAM block inspector has been disabled.");
            case ShowIndex:
                showIndexBins(new File(bamFileName+".bai"));
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

    private void showIndexBins(File bamFile) {
        BAMFileIndexContentInspector inspector = new BAMFileIndexContentInspector(bamFile);
        inspector.inspect(System.out,null,null);
    }

    private class BAMFileIndexContentInspector /*extends CachingBAMFileIndex*/ {
        public BAMFileIndexContentInspector(File bamFileIndex) {
//            super(bamFileIndex);    
        }

        public void inspect(PrintStream outputStream, Integer startPosition, Integer stopPosition) {
            /*
            outputStream.printf("Number of reference sequences: %d%n", this.referenceToBins.size());
            for(int referenceSequence: referenceToBins.keySet()) {
                Bin[] bins = referenceToBins.get(referenceSequence);
                outputStream.printf("Reference sequence: %d%n",referenceSequence);
                outputStream.printf("number of bins: %d%n",bins.length);
                for(Bin bin: bins) {
                    List<Chunk> chunks = binToChunks.get(bin);
                    outputStream.printf("\tBin: %d, number of chunks: %d%n", bin.binNumber, chunks.size());
                    for(Chunk chunk: chunks)
                        outputStream.printf("\t\tChunk: %s%n", chunk);
                }
                LinearIndex linearIndex = referenceToLinearIndices.get(referenceSequence);
                outputStream.printf("\t\tIndex entries: %d", linearIndex.indexEntries.length);
                for(long indexEntry: linearIndex.indexEntries)
                    outputStream.printf("%d,",indexEntry);
                outputStream.printf("%n");
            }
            */
        }
    }
}
