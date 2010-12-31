/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.qc;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.util.ParsingUtils;
import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.vcf.VCFCodec;
import org.broad.tribble.readers.AsciiLineReader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.SimpleTimer;

import java.io.PrintStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.*;

/**
 * Emits specific fields as dictated by the user from one or more VCF files.
 */
@Requires(value={})
public class ProfileRodSystem extends RodWalker<Integer, Integer> {
    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(fullName="nIterations", shortName="N", doc="Number of raw reading iterations to perform", required=false)
    int nIterations = 1;

    @Argument(fullName="maxRecords", shortName="M", doc="Max. number of records to process", required=false)
    int MAX_RECORDS = -1;

    SimpleTimer timer = new SimpleTimer("myTimer");

    public void initialize() {
        File rodFile = getRodFile();

        out.printf("# walltime is in seconds%n");
        out.printf("# file is %s%n", rodFile);
        out.printf("# file size is %d bytes%n", rodFile.length());
        out.printf("operation\titeration\twalltime%n");
        for ( int i = 0; i < nIterations; i++ ) {
            out.printf("read.bytes\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.BY_BYTE));
            out.printf("read.line\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.BY_LINE));
            out.printf("line.and.parts\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.BY_PARTS));
            out.printf("decode.loc\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.DECODE_LOC));
            out.printf("full.decode\t%d\t%.2f%n", i, readFile(rodFile, ReadMode.DECODE));
        }

        timer.start(); // start up timer for map itself
    }

    private enum ReadMode { BY_BYTE, BY_LINE, BY_PARTS, DECODE_LOC, DECODE };

    private final double readFile(File f, ReadMode mode) {
        timer.start();

        try {
            byte[] data = new byte[100000];
            FileInputStream s = new FileInputStream(f);

            if ( mode == ReadMode.BY_BYTE ) {
                while (true) {
                    if ( s.read(data) == -1 )
                        break;
                }
            } else {
                int counter = 0;
                VCFCodec codec = new VCFCodec();
                String[] parts = new String[100000];
                AsciiLineReader lineReader = new AsciiLineReader(s);

                if ( mode == ReadMode.DECODE_LOC || mode == ReadMode.DECODE )
                    codec.readHeader(lineReader);

                while (counter++ < MAX_RECORDS || MAX_RECORDS == -1) {
                    String line = lineReader.readLine();
                    if ( line == null )
                        break;
                    else if ( mode == ReadMode.BY_PARTS ) {
                        ParsingUtils.split(line, parts, VCFConstants.FIELD_SEPARATOR_CHAR);
                    }
                    else if ( mode == ReadMode.DECODE_LOC ) {
                        codec.decodeLoc(line);
                    }
                    else if ( mode == ReadMode.DECODE ) {
                        processOneVC((VariantContext)codec.decode(line));
                    }
                }
            }
        } catch ( Exception e ) {
            throw new RuntimeException(e);
        }

        return timer.getElapsedTime();
    }

    private File getRodFile() {
        List<ReferenceOrderedDataSource> rods = this.getToolkit().getRodDataSources();
        ReferenceOrderedDataSource rod = rods.get(0);
        return rod.getFile();
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        VariantContext vc = tracker.getVariantContext(ref, "rod", context.getLocation());
        processOneVC(vc);

        return 0;
    }

    private static final void processOneVC(VariantContext vc) {
        vc.getNSamples(); // force us to parse the samples
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        out.printf("gatk.traversal\t%d\t%.2f%n", 0, timer.getElapsedTime());
    }
}