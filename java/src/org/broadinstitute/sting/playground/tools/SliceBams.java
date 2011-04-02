/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.sting.playground.tools;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.picard.util.Log;
import net.sf.samtools.*;
import net.sf.samtools.util.BlockCompressedOutputStream;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.SimplifyingSAMFileWriter;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

/**
 * Reads a list of BAM files and slices all of them into a single merged BAM file
 * containing reads in overlapping chr:start-stop interval.
 *
 * @author Mark DePristo
 */
public class SliceBams extends CommandLineProgram {
    private static final Log log = Log.getInstance(SliceBams.class);

    // Usage and parameters
    @Usage
    public String USAGE = "Merges multiple SAM/BAM files into one BAM overlapping chr:start-stop interval .\n";

    @Option(shortName="I", doc="List of input BAM files")
    public File INPUT_LIST;

    @Option(shortName="O", doc="SAM or BAM file to write merged result to")
    public File OUTPUT;

    @Option(shortName="L", doc="Location to include")
    public String SLICE;

    private static final int PROGRESS_INTERVAL = 1000000;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new SliceBams().instanceMain(argv));
    }

    private List<File> parseInputFiles(File list) {
        try {
            final List<File> files = new ArrayList<File>();
            for (String fileName : new XReadLines(list).readLines() ) {
                files.add(new File(fileName));
            }
            return files;
        } catch ( FileNotFoundException e ) {
            throw new PicardException("Couldn't read input list", e);
        }
    }

    /**
     * Walk over the input files, reading the headers, and finally prepare the output
     * BAM containing a merge of all of the headers.
     *
     * @param inputBAMs
     * @return
     */
    private SAMFileWriter createOutputBAM(List<File> inputBAMs) {
        SAMFileHeader header = null;

        log.info("Reading headers");
        int fileCounter = 1;
        for (final File inFile : inputBAMs) {
            IoUtil.assertFileIsReadable(inFile);
            final SAMFileReader inReader = new SAMFileReader(inFile, null); // null because we don't want it to look for the index
            final SAMFileHeader inHeader = inReader.getFileHeader();
            log.info("  Reading header from file " + inFile + " " + fileCounter++ + " of " + inputBAMs.size());

            if (header == null) {
                header = inHeader;
            }
            else {
                for ( SAMReadGroupRecord rg : inHeader.getReadGroups() )
                    header.addReadGroup(rg);
            }

            inReader.close();
        }

        SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);
        return new SimplifyingSAMFileWriter(out);
    }

    /** Combines multiple SAM/BAM files into one. */
    @Override
	protected int doWork() {
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);

        // Open the files for reading and writing
        List<File> inputBAMs = parseInputFiles(INPUT_LIST);
        IoUtil.assertFileIsWritable(OUTPUT);
        final SAMFileWriter out = createOutputBAM(inputBAMs);
        GenomeLocParser glParser = new GenomeLocParser(out.getFileHeader().getSequenceDictionary());
        GenomeLoc loc = glParser.parseGenomeLoc(SLICE);

        log.info("Reading BAM records");
        long numRecords = 1;
        int fileCounter = 1;
        for (final File inFile : inputBAMs) {
            IoUtil.assertFileIsReadable(inFile);
            log.info("  Reading file " + inFile + " " + fileCounter++ + " of " + inputBAMs.size());
            final SAMFileReader reader = new SAMFileReader(inFile);
            SAMRecordIterator iterator = reader.queryOverlapping(loc.getContig(), loc.getStart(), loc.getStop());

            while ( iterator.hasNext() ) {
                final SAMRecord record = iterator.next();
                out.addAlignment(record);
                if (numRecords % PROGRESS_INTERVAL == 0) {
                    log.info(numRecords + " records read.");
                }
            }

            reader.close();
        }

        log.info("Finished reading inputs.");
        log.info("Sorting final output file.");
        out.close();
        return 0;
    }
}
