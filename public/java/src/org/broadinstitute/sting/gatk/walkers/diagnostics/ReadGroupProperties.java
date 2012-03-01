/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.diagnostics;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.Median;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

/**
 * Emits a GATKReport containing read group, sample, library, platform, center, median insert size and
 * median read length for each read group in every BAM file.
 *
 * Note that this walker stops when all read groups have been observed at least a few thousand times so that
 * the median statistics are well determined.  It is safe to run it WG and it'll finish in an appropriate
 * timeframe.
 *
 * <h2>Input</h2>
 *  <p>
 *      Any number of BAM files
 *  </p>
 *
 * <h2>Output</h2>
 *  <p>
 *      GATKReport containing read group, sample, library, platform, center, median insert size and median read length.
 *
 *      For example, running this tool on the NA12878 data sets:
 *
 *      <pre>
 *      ##:GATKReport.v0.2 ReadGroupProperties : Table of read group properties
 *      readgroup  sample   library       platform  center  median.read.length  median.insert.size
 *      20FUK.1    NA12878  Solexa-18483  illumina  BI                     101                 387
 *      20FUK.2    NA12878  Solexa-18484  illumina  BI                     101                 415
 *      20FUK.3    NA12878  Solexa-18483  illumina  BI                     101                 388
 *      20FUK.4    NA12878  Solexa-18484  illumina  BI                     101                 415
 *      20FUK.5    NA12878  Solexa-18483  illumina  BI                     101                 387
 *      20FUK.6    NA12878  Solexa-18484  illumina  BI                     101                 415
 *      20FUK.7    NA12878  Solexa-18483  illumina  BI                     101                 388
 *      20FUK.8    NA12878  Solexa-18484  illumina  BI                     101                 415
 *      20GAV.1    NA12878  Solexa-18483  illumina  BI                     101                 388
 *      20GAV.2    NA12878  Solexa-18484  illumina  BI                     101                 415
 *      20GAV.3    NA12878  Solexa-18483  illumina  BI                     101                 388
 *      20GAV.4    NA12878  Solexa-18484  illumina  BI                     101                 416
 *      20GAV.5    NA12878  Solexa-18483  illumina  BI                     101                 388
 *      20GAV.6    NA12878  Solexa-18484  illumina  BI                     101                 415
 *      20GAV.7    NA12878  Solexa-18483  illumina  BI                     101                 387
 *      20GAV.8    NA12878  Solexa-18484  illumina  BI                     101                 414
 *      </pre>
 *  </p>
 *
 * <h2>Examples</h2>
 *  <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T ReadGroupProperties
 *      -I example1.bam -I example2.bam etc
 *      -R reference.fasta
 *      -o example.gatkreport.txt
 *  </pre>
 *
 * @author Mark DePristo
 */



public class ReadGroupProperties extends ReadWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @Argument(shortName="maxElementsForMedian", doc="Calculate median from the first maxElementsForMedian values observed", required=false)
    public int MAX_VALUES_FOR_MEDIAN = 10000;

    private final static String TABLE_NAME = "ReadGroupProperties";
    private final Map<String, Median<Integer>> readLengths = new HashMap<String, Median<Integer>>();
    private final Map<String, Median<Integer>> insertSizes = new HashMap<String, Median<Integer>>();

    @Override
    public void initialize() {
        for ( final SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
            readLengths.put(rg.getId(), new Median<Integer>(MAX_VALUES_FOR_MEDIAN));
            insertSizes.put(rg.getId(), new Median<Integer>(MAX_VALUES_FOR_MEDIAN));
        }
    }

    @Override
    public boolean filter(ReferenceContext ref, GATKSAMRecord read) {
        return ! (read.getReadFailsVendorQualityCheckFlag() || read.getReadUnmappedFlag());
    }

    @Override
    public boolean isDone() {
        // TODO -- this is far too slow!
        return ! (anyMedianNeedsData(readLengths) || anyMedianNeedsData(insertSizes));
    }

    private final boolean anyMedianNeedsData(Map<String, Median<Integer>> medianMap) {
        for ( Median<Integer> median : medianMap.values() ) {
            if ( ! median.isFull() )
                return true;
        }

        return false;
    }

    private final void updateMedian(final Median<Integer> median, final int value) {
        median.add(value);
    }

    @Override
    public Integer map(ReferenceContext referenceContext, GATKSAMRecord read, ReadMetaDataTracker readMetaDataTracker) {
        final String rg = read.getReadGroup().getId();

        updateMedian(readLengths.get(rg), read.getReadLength());
        if ( read.getReadPairedFlag() && read.getInferredInsertSize() != 0) {
            //logger.info(rg + " => " + Math.abs(read.getInferredInsertSize()));
            updateMedian(insertSizes.get(rg), Math.abs(read.getInferredInsertSize()));
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer integer, Integer integer1) {
        return null;
    }

    @Override
    public void onTraversalDone(Integer sum) {
        final GATKReport report = new GATKReport();
        report.addTable(TABLE_NAME, "Table of read group properties");
        GATKReportTable table = report.getTable(TABLE_NAME);

        table.addPrimaryKey("readgroup");
        //* Emits a GATKReport containing read group, sample, library, platform, center, median insert size and
        //* median read length for each read group in every BAM file.
        table.addColumn("sample", "NA");
        table.addColumn("library", "NA");
        table.addColumn("platform", "NA");
        table.addColumn("center", "NA");
        table.addColumn("median.read.length", Integer.valueOf(0));
        table.addColumn("median.insert.size", Integer.valueOf(0));

        for ( final SAMReadGroupRecord rg : getToolkit().getSAMFileHeader().getReadGroups() ) {
            final String rgID = rg.getId();
            table.set(rgID, "sample", rg.getSample());
            table.set(rgID, "library", rg.getLibrary());
            table.set(rgID, "platform", rg.getPlatform());
            table.set(rgID, "center", rg.getSequencingCenter());
            table.set(rgID, "median.read.length", readLengths.get(rgID).getMedian(0));
            table.set(rgID, "median.insert.size", insertSizes.get(rgID).getMedian(0));
        }

        report.print(out);
    }
}
