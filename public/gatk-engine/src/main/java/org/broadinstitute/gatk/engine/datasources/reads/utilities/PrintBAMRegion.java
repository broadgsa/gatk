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

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 25, 2011
 * Time: 3:25:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class PrintBAMRegion extends CommandLineProgram {
    @Argument(fullName="input",shortName="I",doc="Input file to process",required=true)
    private File input = null;

    @Argument(fullName="region",shortName="R",doc="BAM region to process, in chunk format (mmmm:nn-xxxx:yy)",required=true)
    private String region;

    private static final long MIN_BLOCK_SIZE = 0;
    private static final long MAX_BLOCK_SIZE = (long)Math.pow(2,48)-1;
    private static final int MIN_OFFSET_SIZE = 0;
    private static final int MAX_OFFSET_SIZE = (int)Math.pow(2,16)-1;

    public int execute() {
        SAMFileReader reader = new SAMFileReader(input);
        reader.setValidationStringency(ValidationStringency.SILENT);

        Pattern regionPattern = Pattern.compile("(\\d+):(\\d+)-(\\d+):(\\d+)");
        Matcher matcher = regionPattern.matcher(region);
        if(!matcher.matches())
            throw new UserException("BAM region to process must be in chunk format (mmmm:nn-xxxx:yy)");

        long firstBlock = Long.parseLong(matcher.group(1));
        int firstOffset = Integer.parseInt(matcher.group(2));
        long lastBlock = Long.parseLong(matcher.group(3));
        int lastOffset = Integer.parseInt(matcher.group(4));

        if(firstBlock < MIN_BLOCK_SIZE || firstBlock > MAX_BLOCK_SIZE)
            throw new UserException(String.format("First block is invalid; must be between %d and %d; actually is %d",MIN_BLOCK_SIZE,MAX_BLOCK_SIZE,firstBlock));
        if(lastBlock < MIN_BLOCK_SIZE || lastBlock > MAX_BLOCK_SIZE)
            throw new UserException(String.format("Last block is invalid; must be between %d and %d; actually is %d",MIN_BLOCK_SIZE,MAX_BLOCK_SIZE,lastBlock));
        if(firstOffset < MIN_OFFSET_SIZE || firstOffset > MAX_OFFSET_SIZE)
            throw new UserException(String.format("First offset is invalid; must be between %d and %d; actually is %d",MIN_OFFSET_SIZE,MAX_OFFSET_SIZE,firstOffset));
        if(lastOffset < MIN_OFFSET_SIZE || lastOffset > MAX_OFFSET_SIZE)
            throw new UserException(String.format("Last offset is invalid; must be between %d and %d; actually is %d",MIN_OFFSET_SIZE,MAX_OFFSET_SIZE,lastOffset));

        GATKChunk chunk = new GATKChunk(firstBlock<<16 | firstOffset,lastBlock<<16 | lastOffset);
        GATKBAMFileSpan fileSpan = new GATKBAMFileSpan(chunk);

        SAMRecordIterator iterator = reader.iterator(fileSpan);
        long readCount = 0;
        while(iterator.hasNext()) {
            System.out.printf("%s%n",iterator.next().format());
            readCount++;
        }
        System.out.printf("%d reads shown.",readCount);

        iterator.close();
        reader.close();

        return 0;
    }


    /**
     * Required main method implementation.
     * @param argv Command-line argument text.
     * @throws Exception on error.
     */
    public static void main(String[] argv) throws Exception {
        try {
            PrintBAMRegion instance = new PrintBAMRegion();
            start(instance, argv);
            System.exit(0);
        }
        catch(Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }
}
