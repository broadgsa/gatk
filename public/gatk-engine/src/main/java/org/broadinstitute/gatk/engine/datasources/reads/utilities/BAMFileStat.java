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

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.instrumentation.Sizeof;

import java.io.File;
import java.lang.reflect.Field;
import java.util.List;
import java.util.Map;

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

    @Argument(doc="The range to inspect.",required=false)
    private String range;

    public int execute() {
        switch(command) {
            case ShowBlocks:
                throw new ReviewedGATKException("The BAM block inspector has been disabled.");
            case ShowIndex:
                showIndexBins(new File(bamFileName),range);
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

    private void showIndexBins(File bamFile,String contigName) {
        SAMFileReader reader;
        BAMIndex index;

        reader = new SAMFileReader(bamFile);
        reader.setValidationStringency(ValidationStringency.SILENT);
        reader.enableIndexCaching(true);
        index = reader.getIndex();

        reader.queryOverlapping(contigName,1,reader.getFileHeader().getSequence(contigName).getSequenceLength()).close();

        int numBins = 0;
        int numChunks = 0;
        int numLinearIndexEntries = 0;

        try {
            Field[] fields = index.getClass().getDeclaredFields();
            for(Field field: fields) {
                if(field.getName().equals("mLastReferenceRetrieved")) {
                    field.setAccessible(true);
                    Integer lastReferenceRetrieved = (Integer)field.get(index);
                    System.out.printf("Last reference retrieved: %d%n", lastReferenceRetrieved);
                }

                if(field.getName().equals("mQueriesByReference")) {
                    field.setAccessible(true);
                    Map<Integer,Object> cachedQueries = (Map<Integer,Object>)field.get(index);

                    for(Object bamIndexContent: cachedQueries.values()) {
                        List<Object> bins = null;
                        Map<Object,Object> binToChunkMap = null;
                        Object linearIndex = null;

                        Field[] indexContentFields = bamIndexContent.getClass().getDeclaredFields();
                        for(Field indexContentField: indexContentFields) {
                            if(indexContentField.getName().equals("mReferenceSequence")) {
                                indexContentField.setAccessible(true);
                                System.out.printf("Reference sequence: %d%n", indexContentField.getInt(bamIndexContent));
                            }

                            if(indexContentField.getName().equals("mBins")) {
                                indexContentField.setAccessible(true);
                                bins = (List<Object>)indexContentField.get(bamIndexContent);
                            }

                            if(indexContentField.getName().equals("mBinToChunks")) {
                                indexContentField.setAccessible(true);
                                binToChunkMap = (Map<Object,Object>)indexContentField.get(bamIndexContent);
                            }

                            if(indexContentField.getName().equals("mLinearIndex")) {
                                indexContentField.setAccessible(true);
                                linearIndex = indexContentField.get(bamIndexContent);
                            }
                        }

                        numBins = bins.size();
                        for(Object bin: bins) {
                            int binNumber;

                            Field[] binFields = bin.getClass().getDeclaredFields();
                            for(Field binField: binFields) {
                                if(binField.getName().equals("binNumber")) {
                                    binField.setAccessible(true);
                                    binNumber = binField.getInt(bin);
                                    List<Object> chunks = (List<Object>)binToChunkMap.get(bin);
                                    System.out.printf("\tBin: %d, number of chunks: %d%n",binNumber,chunks.size());
                                    for(Object chunk: chunks)
                                        System.out.printf("\t\tChunk: %s%n",chunk);
                                    numChunks += chunks.size();
                                }
                            }
                        }

                        Field[] linearIndexFields = linearIndex.getClass().getDeclaredFields();
                        for(Field linearIndexField: linearIndexFields) {
                            if(linearIndexField.getName().equals("mIndexEntries")) {
                                linearIndexField.setAccessible(true);
                                long[] linearIndexEntries = (long[])linearIndexField.get(linearIndex);
                                System.out.printf("\t\tIndex entries: %d", linearIndexEntries.length);
                                for(long indexEntry: linearIndexEntries)
                                    System.out.printf("%d,",indexEntry);
                                System.out.printf("%n");
                                numLinearIndexEntries = linearIndexEntries.length;
                            }
                        }
                    }
                }
            }
        }
        catch(IllegalAccessException ex) {
            throw new ReviewedGATKException("Unable to examine cached index",ex);
        }

        System.out.printf("%nOverall: %d bins, %d chunks, %d linear index entries",numBins,numChunks,numLinearIndexEntries);
        if(Sizeof.isEnabled())
            System.out.printf(", total index size in bytes: %d",Sizeof.getObjectGraphSize(index));
        System.out.println();

        reader.close();
    }
}
