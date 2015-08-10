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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.datasources.reads.FilePointer;
import org.broadinstitute.gatk.engine.datasources.reads.IntervalSharder;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import org.broadinstitute.gatk.utils.text.ListFileUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

/**
 * Traverses a region in a dataset looking for outliers.
 */
public class FindLargeShards extends CommandLineProgram {
    private static Logger logger = Logger.getLogger(FindLargeShards.class);

    @Input(fullName = "input_file", shortName = "I", doc = "SAM or BAM file(s)", required = false)
    public List<String> samFiles = new ArrayList<String>();

    @Input(fullName = "reference_sequence", shortName = "R", doc = "Reference sequence file", required = false)
    public File referenceFile = null;

    @Input(fullName = "intervals", shortName = "L", doc = "A list of genomic intervals over which to operate. Can be explicitly specified on the command line or in a file.",required=false)
    public List<String> intervals = null;

    @Output(required=false)
    public PrintStream out = System.out;

    /**
     * The square of the sum of all uncompressed data.  Based on the BAM spec, the size of this could be
     * up to (2^64)^2.
     */
    private BigInteger sumOfSquares = BigInteger.valueOf(0);

    /**
     * The running sum of all uncompressed data.  Based on the BAM spec, the BAM must be less than Long.MAX_LONG
     * when compressed -- in other words, the sum of the sizes of all BGZF blocks must be < 2^64.
     */
    private BigInteger sum = BigInteger.valueOf(0);

    /**
     * The number of shards viewed.
     */
    private long numberOfShards;


    @Override
    public int execute() throws IOException {
        // initialize reference
        IndexedFastaSequenceFile refReader = new IndexedFastaSequenceFile(referenceFile);
        GenomeLocParser genomeLocParser = new GenomeLocParser(refReader);        

        // initialize reads
        List<SAMReaderID> bamReaders = ListFileUtils.unpackBAMFileList(samFiles,parser);
        SAMDataSource dataSource = new SAMDataSource(referenceFile, bamReaders, new ThreadAllocation(), null, genomeLocParser);

        // intervals
        final GenomeLocSortedSet intervalSortedSet;
        if ( intervals != null )
            intervalSortedSet = IntervalUtils.sortAndMergeIntervals(genomeLocParser, IntervalUtils.parseIntervalArguments(genomeLocParser, intervals), IntervalMergingRule.ALL);
        else
            intervalSortedSet = GenomeLocSortedSet.createSetFromSequenceDictionary(refReader.getSequenceDictionary());

        logger.info(String.format("PROGRESS: Calculating mean and variance: Contig\tRegion.Start\tRegion.Stop\tSize"));        

        IntervalSharder sharder = IntervalSharder.shardOverIntervals(dataSource,intervalSortedSet,IntervalMergingRule.ALL);
        while(sharder.hasNext()) {
            FilePointer filePointer = sharder.next();

            // Size of the file pointer.
            final long size = filePointer.size();            

            BigInteger bigSize = BigInteger.valueOf(size);
            sumOfSquares = sumOfSquares.add(bigSize.pow(2));
            sum = sum.add(bigSize);
            numberOfShards++;

            if(numberOfShards % 1000 == 0) {
                GenomeLoc boundingRegion = getBoundingRegion(filePointer,genomeLocParser);
                logger.info(String.format("PROGRESS: Calculating mean and variance: %s\t%d\t%d\t%d",boundingRegion.getContig(),boundingRegion.getStart(),boundingRegion.getStop(),size));
            }

        }

        // Print out the stddev: (sum(x^2) - (1/N)*sum(x)^2)/N
        long mean = sum.divide(BigInteger.valueOf(numberOfShards)).longValue();
        long stddev = (long)(Math.sqrt(sumOfSquares.subtract(sum.pow(2).divide(BigInteger.valueOf(numberOfShards))).divide(BigInteger.valueOf(numberOfShards)).doubleValue()));
        logger.info(String.format("Number of shards: %d; mean uncompressed size = %d; stddev uncompressed size  = %d%n",numberOfShards,mean,stddev));

        // Crank through the shards again, this time reporting on the shards significantly larger than the mean.
        long threshold = mean + stddev*5;
        logger.warn(String.format("PROGRESS: Searching for large shards: Contig\tRegion.Start\tRegion.Stop\tSize"));
        out.printf("Contig\tRegion.Start\tRegion.Stop\tSize%n");

        sharder =  IntervalSharder.shardOverIntervals(dataSource,intervalSortedSet,IntervalMergingRule.ALL);
        while(sharder.hasNext()) {
            FilePointer filePointer = sharder.next();

            // Bounding region.
            GenomeLoc boundingRegion = getBoundingRegion(filePointer,genomeLocParser);

            // Size of the file pointer.
            final long size = filePointer.size();            

            numberOfShards++;

            if(filePointer.size() <= threshold) {
                if(numberOfShards % 1000 == 0) 
                    logger.info(String.format("PROGRESS: Searching for large shards: %s\t%d\t%d\t%d",boundingRegion.getContig(),boundingRegion.getStart(),boundingRegion.getStop(),size));
                continue;
            }

            out.printf("%s\t%d\t%d\t%d%n",boundingRegion.getContig(),boundingRegion.getStart(),boundingRegion.getStop(),size);
        }

        return 0;
    }

    private GenomeLoc getBoundingRegion(final FilePointer filePointer, final GenomeLocParser genomeLocParser) {
        List<GenomeLoc> regions = filePointer.getLocations();

        // The region contained by this FilePointer.
        final String contig = regions.get(0).getContig();
        final int start = regions.get(0).getStart();
        final int stop = regions.get(regions.size()-1).getStop();

        return genomeLocParser.createGenomeLoc(contig,start,stop);
    }

    /**
     * Required main method implementation.
     * @param argv Command-line argument text.
     * @throws Exception on error.
     */
    public static void main(String[] argv) throws Exception {
        int returnCode = 0;
        try {
            FindLargeShards instance = new FindLargeShards();
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
